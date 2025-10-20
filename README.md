# MockDE

**MockDE** is a lightweight FastAPI backend that simulates differential gene expression (DE) analysis using mock RNA-seq data. It provides a `/diffexp` endpoint for testing DE workflows without requiring real TCGA datasets.

---

## 🚀 Quick Start

### Requirements

- **Python 3.11**
- Recommended: virtual environment

---

### 1. Clone and enter project

```bash
git clone https://github.com/aminbayat3/DifferentialExpressionMockApi.git
cd DifferentialExpressionMockApi
```

---

### 2. Create and activate virtual environment

#### 🪟 Windows (PowerShell)

```bash
python -m venv venv
venv\Scripts\Activate.ps1
```

#### 💻 macOS / Linux

```bash
python3 -m venv venv
source venv/bin/activate
```

---

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

---

### 4. Run the API server

```bash
uvicorn app.main:app --reload
```

Then open your browser at **[http://127.0.0.1:8000/docs](http://127.0.0.1:8000/docs)** to test the `/diffexp` endpoint.

---
