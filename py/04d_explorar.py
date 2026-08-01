"""
04d — Diagnóstico del esquema de emisiones en DiscoData.

No descarga datos: imprime las columnas y una fila de ejemplo de las tablas
implicadas, para poder diseñar la unión con precisión en lugar de adivinar.

Las vistas Tableau_ resultaron incompletas, así que hay que reconstruir la
cadena normalizada:
    PollutantRelease -> ProductionFacilityReport -> ProductionFacility
                                                 -> ProductionSite
Cada eslabón aporta una pieza: la masa emitida, el año y el país, y el
identificador de sitio que enlaza con nuestro patrón de puntos.
"""

from __future__ import annotations

import json
from urllib.parse import quote

import requests

from common import RAW, log

URL_SQL = "https://discodata.eea.europa.eu/sql"
BASE = "[IED].[v1r2]"

TABLAS = [
    "EmissionsToAir",
    "Tableau_Browse6_PollutantReleases",
    "ProductionFacilityReport",
    "ProductionFacility",
    "ProductionSite",
    "PollutantRelease",
]


def consultar(sql: str, n: int = 1) -> list[dict]:
    url = f"{URL_SQL}?query={quote(sql)}&p=1&nrOfHits={n}"
    r = requests.get(url, timeout=180)
    r.raise_for_status()
    j = r.json()
    if "errors" in j:
        raise RuntimeError(j["errors"][0].get("error", j["errors"]))
    return j.get("results", [])


def recortar(v, n: int = 44) -> str:
    s = str(v)
    return s if len(s) <= n else s[: n - 1] + "…"


def main() -> None:
    esquema = {}
    for t in TABLAS:
        print("\n" + "=" * 74)
        print(f"  {t}")
        print("=" * 74)
        try:
            filas = consultar(f"SELECT TOP 1 * FROM {BASE}.[{t}]")
        except Exception as exc:  # noqa: BLE001
            print(f"  error: {str(exc)[:90]}")
            continue
        if not filas:
            print("  (sin filas)")
            continue
        fila = filas[0]
        esquema[t] = list(fila.keys())
        for k, v in fila.items():
            print(f"    {k:<34} = {recortar(v)}")

    (RAW / "ied_esquema_emisiones.json").write_text(
        json.dumps(esquema, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # Un dato clave: cómo se escribe el NOx en la columna de contaminante.
    print("\n" + "=" * 74)
    print("  VALORES DE 'pollutant' QUE CONTIENEN NITRÓGENO")
    print("=" * 74)
    for t in ("EmissionsToAir", "PollutantRelease"):
        if t not in esquema:
            continue
        col = next((c for c in esquema[t] if c.lower() == "pollutant"), None)
        if not col:
            continue
        try:
            d = consultar(
                f"SELECT DISTINCT TOP 80 {col} AS v FROM {BASE}.[{t}]", n=80
            )
            vals = sorted({str(x["v"]) for x in d if x.get("v")})
            print(f"\n  {t}  ({len(vals)} valores distintos):")
            for v in vals:
                marca = "  <---" if ("itrogen" in v or "NOX" in v.upper()) else ""
                print(f"      {v}{marca}")
        except Exception as exc:  # noqa: BLE001
            print(f"  {t}: {str(exc)[:80]}")

    log("\nEsquema guardado en data/raw/ied_esquema_emisiones.json")


if __name__ == "__main__":
    main()
