#!/bin/bash

# Define paths
WORKSPACE_DIR="/home/sriman/Documents/Thesis/workspace"
RAD_DIR="${WORKSPACE_DIR}/programs/rad"
SRC_FILE="${RAD_DIR}/src/something.f90"
INPUT_FILE="${RAD_DIR}/Energy_states_2s_2p.txt"
EXEC_NAME="${RAD_DIR}/bin/something_exec"
INCLUDE_DIR="${WORKSPACE_DIR}/include"
MOD_DIR="${INCLUDE_DIR}/mod"
LIB_DIR="${INCLUDE_DIR}/lib"

# --- Compilation ---
echo "Compiling ${SRC_FILE}..."
# Adjust libraries (-lbspqed, -lslatec etc.) as needed based on your actual dependencies
gfortran -I"${MOD_DIR}" "${SRC_FILE}" -o "${EXEC_NAME}" -L"${LIB_DIR}" -lbspqed -lslatec -static

# Check if compilation was successful
if [ $? -ne 0 ]; then
  echo "Compilation failed. Exiting."
  exit 1
fi
echo "Compilation successful. Executable: ${EXEC_NAME}"

# --- Execution Loop ---
echo "Starting batch processing of ${INPUT_FILE}..."

# Read the input file line by line, skipping the header
tail -n +2 "${INPUT_FILE}" | while IFS= read -r line || [[ -n "$line" ]]; do
  # Skip empty lines potentially introduced by tail or read
  if [[ -z "$line" ]]; then
      continue
  fi

  # Extract values using awk for robustness with spacing
  Z=$(echo "$line" | awk '{print $1}')
  N=$(echo "$line" | awk '{print $2}')
  R_fm=$(echo "$line" | awk '{print $3}')
  E_initial=$(echo "$line" | awk '{print $4}')
  E_intermediate=$(echo "$line" | awk '{print $5}')

  # Print progress
  echo "Processing Z=${Z}, N=${N}, R_fm=${R_fm}..."

  # Prepare the input string for the Fortran program
  # Order=2, Z, N, R_fm, InitialState='2s', IntermediateState='2p', E_initial, E_intermediate
  input_string=$(printf "2\n%s\n%s\n%s\n2s\n2p\n%s\n%s\n" \
                 "${Z}" "${N}" "${R_fm}" "${E_initial}" "${E_intermediate}")

  # Run the executable, piping the input string
  echo -e "${input_string}" | "${EXEC_NAME}"

  # Optional: Add a small delay if needed, e.g., sleep 0.1
done

echo "Batch processing finished."
echo "Results should be appended to /home/sriman/Documents/Thesis/Z_Matrix_elements.txt"

exit 0