<img width="1235" height="790" alt="Screenshot from 2026-09-09 23-02-58" src="https://github.com/user-attachments/assets/2a81ee85-d572-4018-86b8-ee867ee4ac5e" />

---


### Need Only 4 things installed:

**1. Install curl**
```bash
# Ubuntu/Linux
sudo apt update
sudo apt install curl -y
curl --version
```



**1. Install Java (Nextflow required Java 11+.)**
```bash
# Ubuntu/Linux
sudo apt update
sudo apt install openjdk-17-jdk -y
java -version
```

**1. Assign the appropriate Java version (If Needed)**
```bash
# Ubuntu/Linux
sudo update-alternatives --config java

```



**1. Install Docker**
```bash
# Ubuntu/Linux
sudo apt update
sudo apt install docker.io -y
sudo systemctl start docker
sudo systemctl enable docker
sudo usermod -aG docker $USER
docker --version
# Then restart the PC
```

**2. Install Nextflow**
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/          *(moving system path)
```

---

### Then just 3 commands:

**Step 1 — Clone the project**
```bash
git clone https://github.com/YOUR_USERNAME/Progression_CRC.git
cd Progression_CRC
```

**Step 2 — Build Docker images (one time only, ~15/20 mins)**
```bash
docker build -f Dockerfile.r      -t progression_crc_r:latest      .
docker build -f Dockerfile.python -t progression_crc_python:latest  .
```

**Step 3 — Run the pipeline**
```bash
nextflow run main.nf
```

That's it.

---

### Other run options if needed

```bash
# Resume if pipeline was interrupted
nextflow run main.nf -resume

```

---

### Results

Ensure internet connection during pipeline run. After the pipeline finishes, all outputs are saved automatically:

results/figures/ (All figures)

results/tables/ (All result tables) 

---





