import requests
from bs4 import BeautifulSoup

with open("mibig_ids.txt") as f:
    ids = [x.strip() for x in f if x.strip()]

with open("mibig_exact_products.tsv", "w") as out:
    out.write("MiBIG_ID\tCompound\n")
    for bgc_id in ids:
        url = f"https://mibig.secondarymetabolites.org/go/{bgc_id}"
        compound = "NA"
        try:
            r = requests.get(url, timeout=30)
            if r.ok:
                soup = BeautifulSoup(r.text, "html.parser")
                text = soup.get_text("\n")
                lines = [x.strip() for x in text.splitlines() if x.strip()]
                for i, line in enumerate(lines):
                    if line == "Compounds":
                        vals = []
                        for nxt in lines[i + 1:i + 15]:
                            if nxt in {
                                "Publications",
                                "Notes",
                                "Loci",
                                "Classes",
                                "Biosynthetic class(es)",
                                "NCBI accession"
                            }:
                                break
                            vals.append(nxt)
                        if vals:
                            compound = "; ".join(vals)
                        break
        except Exception:
            compound = "NA"
        out.write(f"{bgc_id}\t{compound}\n")
