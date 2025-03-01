import random

def generate_variant(chrom, pos):
    ref = random.choice("ACGT")
    alt = random.choice([base for base in "ACGT" if base != ref])
    dp = random.randint(100, 1000)
    af = random.uniform(0.001, 0.5)
    ad_alt = int(dp * af)
    ad_ref = dp - ad_alt
    tlod = 0
    confidence = "low"

    if af > 0.2 and dp > 500:
      tlod = random.uniform(50, 1500)  # Higher TLOD for likely real variants
      if af > 0.4:
        confidence = "high"
      else:
        confidence = "medium"
    elif af < 0.1:
        tlod = random.uniform(-2,2) # Lower TLOD for likely not real
    else:
      tlod = random.uniform(-2, 50)


    info = f"AS_SB_TABLE={ad_ref},0|{ad_alt},0;DP={dp};ECNT={random.randint(10, 50)};MBQ=40,40;MFRL=0,0;MMQ=60,60;MPOS={random.randint(5, 100)};POPAF={random.uniform(0.1, 10)};TLOD={tlod}"
    format_fields = f"GT:AD:AF:DP:WC:ML:ND:DC:SBM:PO:MO:SC"
    sample_data = f"0/1:{ad_alt}:{af:.3f}:{dp}:{random.randint(100, 800)}:1:{random.randint(50, 500)}:{random.randint(10, 100)}:0.25:{random.randint(10, 200)}:{random.randint(10, 200)}:{ref}{alt}"

    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\t{format_fields}\t{sample_data}\t{confidence}"


with open("mock_variants.tsv", "w") as f:
    f.write("#tumor_sample=MT5222-061-13_aPCR_Native\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMT5222-061-13_aPCR_Native\tConfidence\n")

    # Initial 4 variants (from your example)
    f.write("chr1\t917495\t.\tC\tT\t.\tPASS\tAS_SB_TABLE=580,0|1,0;DP=590;ECNT=15;MBQ=40,40;MFRL=0,0;MMQ=60,60;MPOS=6;POPAF=7.3;TLOD=-1.056\tGT:AD:AF:DP:WC:ML:ND:DC:SBM:PO:MO:SC\t0/1:1:0.003:534:643:1:421:111:0.249:193:228:GG\tlow\n")
    f.write("chr1\t2556714\t.\tA\tG\t.\tPASS\tAS_SB_TABLE=409,0|377,0;DP=805;ECNT=27;MBQ=40,40;MFRL=0,0;MMQ=60,60;MPOS=51;POPAF=0.274;TLOD=1075.04\tGT:AD:AF:DP:WC:ML:ND:DC:SBM:PO:MO:SC\t0/1:377:0.476:720:409:1:285:62:0.249:134:151:AA\thigh\n")
    f.write("chr1\t2557761\t.\tC\tT\t.\tPASS\tAS_SB_TABLE=565,0|162,0;DP=754;ECNT=38;MBQ=40,40;MFRL=0,0;MMQ=60,60;MPOS=43;POPAF=2.01;TLOD=464.11\tGT:AD:AF:DP:WC:ML:ND:DC:SBM:PO:MO:SC\t0/1:162:0.215:679:171:1:127:22:0.25:65:62:AG\tmedium\n")
    f.write("chr1\t2557872\t.\tG\tA\t.\tPASS\tAS_SB_TABLE=775,0|4,0;DP=791;ECNT=28;MBQ=40,16;MFRL=0,0;MMQ=60,60;MPOS=73;POPAF=7.3;TLOD=-0.9435\tGT:AD:AF:DP:WC:ML:ND:DC:SBM:PO:MO:SC\t0/1:4:0.003:716:5:0:3:1:0.234:2:1:CA\tlow\n")



    # Generate 296 more variants
    for i in range(296):
        pos = 2557872 + i * 1000  # Increment position
        f.write(generate_variant("chr1", pos) + "\n")