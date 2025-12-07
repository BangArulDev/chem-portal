from contextlib import asynccontextmanager
from fastapi import FastAPI, HTTPException, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import Optional

import aiosqlite
import os

# 🟢 Impor RELATIF untuk Agen (menggunakan nama file yang benar)
from .NovelChemicalDiscoveryAgent import NovelChemicalDiscoveryAgent 
from .auth import AuthMiddleware, LoginHandler, RegisterHandler

# ==============================================================================
# 1. ⚙️ KONFIGURASI FASTAPI
# ==============================================================================

@asynccontextmanager
async def lifespan(app: FastAPI):
    app.state.db = await aiosqlite.connect("data/data.db")

    # Initialize user table if it is not exists yet
    await app.state.db.execute("""
        CREATE TABLE IF NOT EXISTS user (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT NOT NULL,
            username TEXT UNIQUE NOT NULL,
            password TEXT NOT NULL
        )
    """)
    await app.state.db.commit()

    yield

    await app.state.db.close()

app = FastAPI(title="Chemical Discovery Agent API", version="1.0.0", lifespan=lifespan)

# Konfigurasi CORS
web_url = os.getenv("FRONTEND_URL")
origins = ["http://127.0.0.1", "http://localhost:8000", "http://localhost:3000", "*"] if web_url is None else [web_url]
print(origins)
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=False,
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["Authorization"],
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

@app.middleware("http")
async def duh(r, c):
    return await AuthMiddleware(r, c)

@app.get("/")
def read_root():
    return {"status": "ok", "message": "Chemical Discovery Agent API is running."}

@app.post("/discover", response_model=DiscoveryResponse)
async def discover_compound(criteria: CriteriaInput):
    if agent.predictor is None or agent.client is None:
        raise HTTPException(status_code=503, detail="Service unavailable.")

    compound_data_raw = agent.predict_and_lookup(criteria.model_dump())
    justification = agent.get_justification(criteria.model_dump(), compound_data_raw)
    
    final_response = {
        "recommended_compound": compound_data_raw["recommended_compound"],
        "justification_ai": justification
    }
    
    return final_response


@app.get("/results/{filename}")
async def get_results_file(filename: str):
    file_path = os.path.join(agent.RESULTS_DIR, filename)
    
    if not os.path.exists(file_path) or not file_path.startswith(agent.RESULTS_DIR):
        raise HTTPException(status_code=404, detail="File not found")
    
    if filename.endswith(".png"):
        media_type = "image/png"
    elif filename.endswith(".mol"):
        media_type = "chemical/x-mol"
    else:
        raise HTTPException(status_code=400, detail="Invalid request")

    return FileResponse(file_path, media_type=media_type)

@app.get("/profile")
async def get_profile(request: Request):
    return {"user": request.state.user}

@app.post("/login")
async def login(req: Request):
    return await LoginHandler(req)

@app.post("/signup")
async def singup(req: Request):
    return await RegisterHandler(req)