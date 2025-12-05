from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import Dict, Any, Optional
import os

# 🟢 Impor RELATIF untuk Agen (menggunakan nama file yang benar)
from .NovelChemicalDiscoveryAgent import NovelChemicalDiscoveryAgent 

# ==============================================================================
# 1. ⚙️ KONFIGURASI FASTAPI
# ==============================================================================

app = FastAPI(title="Chemical Discovery Agent API", version="1.0.0")

# Konfigurasi CORS
origins = ["http://127.0.0.1", "http://localhost:8000", "http://localhost:3000", "*"] 
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Inisialisasi Agen (Memuat model dan Gemini Client saat startup)
agent = NovelChemicalDiscoveryAgent()


# ==============================================================================
# 2. 📚 DEFINISI SKEMA DATA (Pydantic)
# ==============================================================================

class CriteriaInput(BaseModel):
    Solubility: float
    Viscosity_cP: float
    ThermalStability_Score: float
    BoilingPoint_C: float

class CompoundProperties(BaseModel):
    SMILES: str
    IUPAC: str
    InChI: str
    InChIKey: str
    MolWeight: float
    ExactMass: float
    HBondDonors: int
    HBondAcceptors: int
    TPSA: float
    LogP: float
    RotatableBonds: int

class RecommendedCompound(BaseModel):
    SMILES: str
    Properties: CompoundProperties
    Structure_2D_Path: Optional[str]
    Structure_3D_Path: Optional[str]

class DiscoveryResponse(BaseModel):
    recommended_compound: RecommendedCompound
    justification_ai: str


# ==============================================================================
# 3. 🗺️ ENDPOINT API
# ==============================================================================

@app.get("/")
def read_root():
    return {"status": "ok", "message": "Chemical Discovery Agent API is running."}


@app.post("/discover", response_model=DiscoveryResponse)
async def discover_compound(criteria: CriteriaInput):
    if agent.predictor is None or agent.client is None:
        raise HTTPException(status_code=503, detail="Layanan tidak siap.")

    compound_data_raw = agent.predict_and_lookup(criteria.dict())
    justification = agent.get_justification(criteria.dict(), compound_data_raw)
    
    final_response = {
        "recommended_compound": compound_data_raw["recommended_compound"],
        "justification_ai": justification
    }
    
    return final_response


@app.get("/results/{filename}")
async def get_results_file(filename: str):
    file_path = os.path.join(agent.RESULTS_DIR, filename)
    
    if not os.path.exists(file_path) or not file_path.startswith(agent.RESULTS_DIR):
        raise HTTPException(status_code=404, detail="File tidak ditemukan.")
    
    if filename.endswith(".png"):
        media_type = "image/png"
    elif filename.endswith(".mol"):
        media_type = "chemical/x-mol"
    else:
        media_type = "application/octet-stream"

    return FileResponse(file_path, media_type=media_type)