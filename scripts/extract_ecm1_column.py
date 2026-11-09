#!/usr/bin/env python3
"""Extract one published ECM1 column; no invented sediment subdivision.

Usage: python3 scripts/extract_ecm1_column.py crust.zip sediments.zip
Archives: URLs in DATA/reference_models/ecm1_19178_provenance.json.
Only the selected source rows and small tables belong in Git.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile


def extract(path, member, number):
    with ZipFile(path) as archive:
        raw = archive.read(member)
    lines = raw.decode("utf-8-sig").splitlines()
    names = lines[0].split()
    for line in lines[1:]:
        values = line.split()
        if values and int(values[0]) == number:
            return dict(zip(names, values)), names, values, hashlib.sha256(raw).hexdigest()
    raise ValueError(f"ECM1 cell {number} is missing")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("crust_zip", type=Path)
    parser.add_argument("sediments_zip", type=Path)
    args = parser.parse_args()
    out = Path(__file__).resolve().parents[1] / "DATA" / "reference_models"
    out.mkdir(exist_ok=True)
    records = []
    for path, member, short, url in [
        (args.crust_zip, "ECM1/ECM1.txt", "crust", "https://www.earthcrustmodel1.com/_files/archives/b6b6a8_f3827eed2fb947978b38244a7448985b.zip?dn=ECM1.zip"),
        (args.sediments_zip, "ECM1_Sediments/Sed_ECM1.txt", "sediments", "https://www.earthcrustmodel1.com/_files/archives/b6b6a8_5c9fa52fcd5d40c6ad99d43e42187006.zip?dn=ECM1_Sediments.zip"),
    ]:
        row, names, values, source_sha = extract(path, member, 19178)
        assert float(row["Lon"]) == -82.5 and float(row["Lat"]) == 36.5
        (out / f"ecm1_19178_{short}_source.tsv").write_text("\t".join(names) + "\n" + "\t".join(values) + "\n")
        records.append(dict(url=url, archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), member=member, member_sha256=source_sha))
        if short == "crust":
            crust = row
        else:
            sediment = row
    layers = []
    for i in range(1, 4):
        h = float(sediment[f"SedT{i}"]) / 1000  # published sediment thicknesses are metres
        if h > 0:
            layers.append([1, h, *[float(sediment[f"Sed{k}{i}"]) for k in ("VP", "VS", "RHO")]])
    for i in range(1, 4):
        layers.append([2, float(crust[f"TLy{i}"]), *[float(crust[f"{k}{i}"]) for k in ("VP", "VS", "RHO")]])
    # Preserve published Vs/Vp/rho; roundoff correction to last crustal thickness only.
    correction = float(crust["Hc"]) - sum(row[1] for row in layers)
    assert abs(correction) < 0.0002
    layers[-1][1] += correction
    with (out / "ecm1_19178_layers.csv").open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow(["zone", "thickness_km", "Vp_km_s", "Vs_km_s", "density_g_cm3"])
        writer.writerows(layers)
    manifest = dict(cell=19178, longitude=-82.5, latitude=36.5, selection="Published thin-sediment continental Appalachian column; not a station-specific estimate", paper="https://doi.org/10.1016/j.earscirev.2023.104493", sources=records, sediment_thickness_km=sum(row[1] for row in layers if row[0] == 1), moho_km=float(crust["Hc"]), last_crust_thickness_roundoff_correction_km=correction, ecm1_uppermost_mantle={k: float(crust[k]) for k in ("VPN", "VSN", "RHON")}, mantle_choice="AK135 reference mantle replaces ECM1's single uppermost-mantle value; see ak135_mantle.csv and README.md", material_note="Published ECM1 values retained here. The matched synthetic uses these thicknesses/Vs and the explicitly labelled THBI material closure, not these published Vp/rho.")
    (out / "ecm1_19178_provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: manifest[k] for k in ("cell", "sediment_thickness_km", "moho_km")}))


if __name__ == "__main__":
    main()
