import io
from fastapi.testclient import TestClient
from enrichment_network_api.main import app

client = TestClient(app)
csv = "term,genes,adjusted_pvalue\nPathway A,TP53;BRCA1;EGFR,0.001\nPathway B,TP53;MYC,0.02\nPathway C,EGFR;MYC,0.05\n"
r = client.post(
    "/network/build",
    files={"file": ("test.csv", csv, "text/csv")},
    data={"options_json": '{"preset":"enrichr","apply_preset":true}'},
)
assert r.status_code == 200, r.text
data = r.json()
assert data["stats"]["nodes"] > 0
r2 = client.post("/diffusion/bipartite", json={"graph": data["graph"], "top_n": 10})
assert r2.status_code == 200, r2.text
print("Smoke test passed")
