"""Where the external binaries are, resolved rather than defaulted.

WHY THIS FILE EXISTS. Nine scripts carry the line

    default=Path("/home/user/tools/k2pkg/root/share/kraken2-2.17.2/libexec/classify")

which is not a parameter of the measurement; it is the layout of one container, written
down in the place a parameter goes. It is the same fault as a trim of 14 reaching a
library that has no non-template prefix, one level further out: a machine fact as a
default. It survived because every run so far was launched from a directory where that
path happened to exist, so nothing ever exercised the default being wrong. The end-to-end
job on a conda environment exercised it immediately, which is what that job is for.

WHY THE PACKAGE RESOLVES INSTEAD OF THE SCRIPTS BEING EDITED. A stated value always wins,
so a `[tools]` entry is the override and is never checked against anything else. Where
the config states nothing, the package finds the binary and passes it explicitly, so the
scripts' defaults do not fire under `hapax run` at all. Editing the defaults out of the
scripts is PLAN.md step 7, where the scripts stop being subprocesses; doing it now would
change the recorded `script_sha256` of files whose committed manifests describe runs that
cannot be repeated.

WHY `classify` IS NOT JUST A `which` CALL. The cascade calls Kraken2's raw `classify`
rather than the `kraken2` wrapper, because the wrapper re-reads the database for every
invocation and the cascade invokes it per tier. `classify` is not on PATH in any
distribution of Kraken2: the wrapper knows where it is, and the wrapper is a script that
says so in a variable. So the resolution order is PATH, then the wrapper's own answer,
then the two layouts bioconda and the source tarball use - and then a refusal naming
`[tools] classify`, because a cascade that silently skipped its classifier tier would
report a smaller identifiable fraction and no reason for it.
"""
from __future__ import annotations

import re
import shutil
from pathlib import Path

from .config import Config, ConfigError

#: Tool name as `[tools]` and the registry of scripts spell it -> what to look for on
#: PATH. `classify` has no PATH name of its own; see `_find_classify`.
ON_PATH = {"diamond": "diamond", "minimap2": "minimap2", "megahit": "megahit"}


def _find_classify() -> Path | None:
    """Kraken2's raw `classify`, by asking the wrapper where it keeps it."""
    direct = shutil.which("classify")
    if direct:
        return Path(direct)
    wrapper = shutil.which("kraken2")
    if not wrapper:
        return None
    wrapper = Path(wrapper).resolve()
    try:
        text = wrapper.read_text(errors="replace")
    except OSError:
        text = ""
    # The wrapper is a Perl script whose install step substitutes the directory holding
    # the compiled helpers into `my $KRAKEN2_DIR = "...";`. That is the authoritative
    # answer for whatever layout this installation uses.
    for match in re.finditer(r'KRAKEN2_DIR\s*=\s*["\']([^"\']+)["\']', text):
        candidate = Path(match.group(1)) / "classify"
        if candidate.exists():
            return candidate
    prefix = wrapper.parent.parent
    guesses = [prefix / "libexec" / "classify", wrapper.parent / "classify"]
    guesses += sorted(prefix.glob("share/kraken2-*/libexec/classify"))
    for candidate in guesses:
        if candidate.exists():
            return candidate
    return None


def resolve(cfg: Config, name: str) -> Path:
    """Where `name` is on this machine. Refuses rather than returning a guess.

    A value under `[tools]` wins and is not checked for plausibility beyond existing,
    because stating it is the override. Everything else is found here, once, so that one
    run's stages all drive the same binary: two tiers of one cascade run by two builds of
    Kraken2 is an operating point nobody could write down.
    """
    stated = cfg.tools.get(name)
    if stated:
        path = Path(stated).expanduser()
        if not path.exists():
            raise ConfigError(
                f"[tools] names {name} = {stated!r} in {cfg.path} and there is nothing "
                f"there. A stated path is the override, so it is used as it is and never "
                f"searched for.")
        return path
    found = _find_classify() if name == "classify" else (
        shutil.which(ON_PATH[name]) if name in ON_PATH else None)
    if found:
        return Path(found)
    hint = ("`classify` is Kraken2's own helper and is not on PATH in any distribution; "
            "it sits beside the wrapper, in libexec. "
            if name == "classify" else "")
    raise ConfigError(
        f"hapax cannot find {name} on this machine. {hint}Install the pinned versions "
        f"(environment.yml, TOOLS.md) or name it under [tools] in {cfg.path}.")


def resolve_all(cfg: Config, names: tuple[str, ...]) -> dict[str, Path]:
    """Resolve several, collecting every absence into one message.

    One refusal naming three missing tools is a machine to fix once; three runs each
    refusing on the next one is the same fix discovered three times.
    """
    found, missing = {}, []
    for name in names:
        try:
            found[name] = resolve(cfg, name)
        except ConfigError as e:
            missing.append(str(e))
    if missing:
        raise ConfigError("\n".join(missing))
    return found
