#!/bin/bash
# VM Setup Script for P1+P2 Pipeline
# Run this after launching an i3.8xlarge or i3.16xlarge instance

set -e

echo "=== Setting up VM for P1+P2 Pipeline ==="

# Mount NVMe
echo "Mounting NVMe storage..."
sudo mkfs.ext4 -F /dev/nvme0n1 || true
sudo mkdir -p /data
sudo mount /dev/nvme0n1 /data || true
sudo chown ubuntu:ubuntu /data
cd /data

# Install dependencies
echo "Installing dependencies..."
sudo apt update
sudo apt install -y openjdk-17-jdk docker.io awscli pigz

# Install Nextflow
echo "Installing Nextflow..."
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Setup Docker
sudo usermod -aG docker ubuntu

# Clone pipeline
echo "Cloning pipeline..."
git clone https://github.com/MarcusAtMotley/PoCFigureReprocessing.git || true
cd PoCFigureReprocessing
git pull

# Create directories
mkdir -p /data/references /data/fastq /data/results

echo ""
echo "=== Setup Complete ==="
echo ""
echo "Next steps:"
echo "1. Log out and back in (for docker group)"
echo "2. Run: aws configure"
echo "3. Run: ./scripts/vm_download_refs.sh"
echo "4. Run: ./scripts/vm_download_fastq.sh"
echo "5. Run: ./scripts/vm_run_p1_p2.sh"
