#!/bin/bash

# ----------------------------
# Parse arguments
# ----------------------------
while getopts "f:s:" opt; do
  case $opt in
    f) INDEX0_FILE="$OPTARG" ;;
    s) SIZE="$OPTARG" ;;
    *) echo "Usage: $0 -f <file.index0> -s <size>"; exit 1 ;;
  esac
done

if [[ -z "$INDEX0_FILE" || -z "$SIZE" ]]; then
  echo "Usage: $0 -f <file.index0> -s <size>"
  exit 1
fi

if [[ ! -f "$INDEX0_FILE" ]]; then
  echo "Error: input file does not exist: $INDEX0_FILE"
  exit 1
fi

# Derived filenames
K2BP_FILE="${INDEX0_FILE}.k2bp"
K2BP_C_FILE="${K2BP_FILE}.c"

echo "====================================="
echo "Running K2BP pipeline"
echo "Input file: $INDEX0_FILE"
echo "Leaf size:  $SIZE"
echo "====================================="

# ----------------------------
# 1. Build
# ----------------------------
echo "[1] k2bp_build..."
./k2bp_build.x -s "$SIZE" "$INDEX0_FILE" || exit 1

# ----------------------------
# 2. Add info
# ----------------------------
echo "[2] k2bp_addinfo (.k2bp)..."
./k2bp_addinfo.x "$K2BP_FILE" -se -p 0.2 || exit 1

# ----------------------------
# 3. Compress
# ----------------------------
echo "[3] k2bp_compress..."
./k2bp_compress.x "$K2BP_FILE" || exit 1

# ----------------------------
# 4. Add info compressed
# ----------------------------
echo "[4] k2bp_addinfo (.k2bp.c)..."
./k2bp_addinfo.x "$K2BP_C_FILE" -se -p 0.2 || exit 1

# ----------------------------
# 5. Info
# ----------------------------
echo "[5] k2bp_info..."
./k2bp_info.x "$K2BP_FILE"
./k2bp_info.x "$K2BP_C_FILE"

# ----------------------------
# 6. Multiplication benchmark
# ----------------------------
echo "[6] Multiplication benchmark..."

/usr/bin/time -f "Command: %C\nS:%S U:%U E:%e Mem(kb):%M\n" \
./k2bp_mul.x "$K2BP_FILE" "$K2BP_FILE"

/usr/bin/time -f "Command: %C\nS:%S U:%U E:%e Mem(kb):%M\n" \
./k2bp_mul.x "$K2BP_C_FILE" "$K2BP_C_FILE"

echo "====================================="
echo "Pipeline finished."
echo "====================================="
