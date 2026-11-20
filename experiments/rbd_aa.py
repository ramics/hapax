import subprocess,sys
# spike CDS starts at 21563 (1-based) in NC_045512.2
ref={}
name=None
for line in open('/home/user/hapax/data/ww_cryptic/NC_045512.2.fa'):
    if line[0]=='>': name=line.split()[0][1:]; ref[name]=[]
    else: ref[name].append(line.strip())
R=''.join(ref['NC_045512.2'])
CODON={}
bases='TCAG'
aas='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i=0
for a in bases:
    for b in bases:
        for c in bases:
            CODON[a+b+c]=aas[i]; i+=1
def aln(qfile, qname, region):
    out=subprocess.run(['minimap2','-c','--cs','-x','asm20','/home/user/hapax/data/ww_cryptic/NC_045512.2.fa',qfile],capture_output=True,text=True).stdout
    for l in out.splitlines():
        f=l.split('\t')
        if f[0]!=qname: continue
        if not (int(f[7])<region[1] and int(f[8])>region[0]): continue
        cs=[x for x in f if x.startswith('cs:Z:')][0][5:]
        return int(f[7]), cs, f[4]
    return None
def muts(start, cs):
    # walk cs string, yield (ref_pos0, refbase, qbase) and deletions
    import re
    pos=start; res=[]
    for tok in re.findall(r':\d+|\*[a-z]{2}|\+[a-z]+|-[a-z]+', cs):
        if tok[0]==':': pos+=int(tok[1:])
        elif tok[0]=='*': res.append((pos,tok[1].upper(),tok[2].upper())); pos+=1
        elif tok[0]=='-': res.append((pos,tok[1:].upper(),'del')); pos+=len(tok)-1
        elif tok[0]=='+': res.append((pos,'ins',tok[1:].upper()))
    return res
S=21563-1  # 0-based CDS start
def aa_changes(ms):
    # group substitutions by codon
    bycodon={}
    dels=[]
    for p,r,q in ms:
        if q=='del': dels.append((p,r)); continue
        if p<S: continue
        ci=(p-S)//3
        bycodon.setdefault(ci,[]).append((p,r,q))
    out=[]
    for ci in sorted(bycodon):
        cs_start=S+ci*3
        codon=list(R[cs_start:cs_start+3])
        orig=''.join(codon)
        for p,r,q in bycodon[ci]:
            codon[p-cs_start]=q
        new=''.join(codon)
        a1=CODON.get(orig,'?'); a2=CODON.get(new,'?')
        if a1!=a2: out.append(f"{a1}{ci+1}{a2}")
        else: out.append(f"({a1}{ci+1} syn)")
    for p,r in dels:
        out.append(f"del nt{p+1}-{p+len(r)} (aa {(p-S)//3+1})")
    return out
for qf,qn in [('hapax_s85_spike/contigs.fa','contig_3'),('megahit_SRR15202285/final.contigs.fa','k119_87')]:
    a=aln(qf,qn,(22700,23300))
    if not a: print(qf,qn,'no alignment'); continue
    st,cs,strand=a
    ms=muts(st,cs)
    print(f"--- {qf} {qn} (ref start {st+1}, strand {strand})")
    print("   ", ', '.join(aa_changes(ms)))
