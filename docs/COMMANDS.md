# Colab: Commands to Run (command line)

> **JupyterHub users**: Replace `python` with `./py` in all commands below.
> The `./py` wrapper ensures the correct interpreter (`/opt/jupyterhub/bin/python3`).

Run these in order in Colab (e.g. in code cells, or in the Colab **Terminal** if you open one from the menu). Drive must be mounted first.

---

## 1. Mount Drive and go to Drive

```bash
# In a Colab Python cell first (only once per session):
# from google.colab import drive
# drive.mount("/content/drive")

cd /content/drive/MyDrive
```

---

## 2. Clone repo to Drive (so it survives disconnect)

```bash
git clone https://github.com/kevincallan/rna-seq.git
cd rna-seq
```

---

## 3. Install deps and download data (references + FASTQs)

```bash
bash scripts/setup_colab.sh
```

---

## 4. Point results/work to Drive (optional but recommended)

```bash
./py scripts/colab_drive_setup.py
```

---

## 5. Run the pipeline

**Full run (single command):**

```bash
./py scripts/run_pipeline.py --config config/config_colab.yaml run
```

**Lower memory (use 1 thread):**

```bash
./py scripts/run_pipeline.py --config config/config_colab.yaml run --threads 1
```

**Skip BigWig step to save memory/time (steps 0–5, 7–11):**

```bash
./py scripts/run_pipeline.py --config config/config_colab.yaml run --threads 1 --steps 0 1 2 3 4 5 7 8 9 10 11
```

**Chunked run (resume if session dies):**  
Use the same `--run-id` for both chunks. After the first run, note the printed run ID (e.g. `GSE48519_colab_20250208_123456`).

Phase 1 (through mapping):

```bash
./py scripts/run_pipeline.py --config config/config_colab.yaml --run-id MY_RUN_ID run --threads 1 --steps 0 1 2 3 4 5
```

Phase 2 (counts, DE, report) — run in a new session if needed; clone to Drive again, then:

```bash
cd /content/drive/MyDrive/rna-seq
bash scripts/setup_colab.sh
./py scripts/colab_drive_setup.py
./py scripts/run_pipeline.py --config config/config_colab.yaml --run-id MY_RUN_ID run --threads 1 --steps 7 8 9 10 11
```

Replace `MY_RUN_ID` with the same value from phase 1.

---

## Reducing memory use

- Use **`--threads 1`** (or set `threads: 1` in `config/config_colab.yaml`).
- **Skip BigWig** (step 6) if you don’t need tracks: use `--steps 0 1 2 3 4 5 7 8 9 10 11`.
- Use **one trimming method** to cut memory and time:  
  `run --methods none`
- **Close other Colab tabs** to leave more RAM for the runtime.

---

## Keeping the session from dying (idle / disconnect)

- **Colab:** There is no supported way to keep a Colab runtime alive indefinitely when idle. The runtime can be recycled after ~90 minutes of inactivity or when the browser tab is closed. Best approach:
  - **Run in chunks** (phase 1 and phase 2 above) and **save everything to Drive** (clone to Drive + `colab_drive_setup.py`). Then you can reconnect and continue with the same `--run-id`.
  - Keep the tab open and avoid long idle periods if you want one long run.
  - Colab Pro can give longer runtimes; otherwise use a local machine or cloud VM.
- **Local / SSH (real command line):** Use **`tmux`** or **`screen`** so the pipeline keeps running if the SSH connection drops:
  ```bash
  tmux new -s rnaseq
  cd /path/to/rna-seq
  ./py scripts/run_pipeline.py --config config/config.yaml run
  # Detach: Ctrl+B then D. Reattach later: tmux attach -t rnaseq
  ```

---

## Local run (no Colab)

```bash
git clone https://github.com/kevincallan/rna-seq.git
cd rna-seq
# Create conda env and install refs/FASTQs per README, then:
./py scripts/run_pipeline.py --config config/config.yaml run
```

Use `tmux` or `screen` as above if you’re on a remote machine over SSH.
