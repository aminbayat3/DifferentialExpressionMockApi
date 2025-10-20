from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.controllers.diffexp_controller import router as diffexp_router

app = FastAPI(title="TCGA DE Backend (Mockable)", version="0.7.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(diffexp_router)