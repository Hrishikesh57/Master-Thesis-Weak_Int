#!/bin/bash

# Path to the executable (!! CORRECTED based on Makefile !!)
EXECUTABLE="/home/sriman/Documents/Thesis/workspace/programs/rad/bin/practise"

# Path to the input file with radii
RADII_FILE="/home/sriman/Documents/Thesis/workspace/programs/rad/radii_ref.txt"

# Directory where the Fortran program is run (and where output files will be created)
RUN_DIR="/home/sriman/Documents/Thesis/workspace/programs/rad"

# Check if the executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    echo "Please check the EXECUTABLE path in the script and ensure the program is compiled."
    echo "You might need to run 'make practise' in the $RUN_DIR directory first."
    exit 1
fi

# Check if the radii file exists
if [ ! -f "$RADII_FILE" ]; then
    echo "Error: Radii file not found at $RADII_FILE"
    exit 1
fi

# Change to the run directory
cd "$RUN_DIR" || { echo "Error: Could not change to directory $RUN_DIR"; exit 1; }

echo "Starting runs..."

# Read the radii file line by line, skipping the header
tail -n +2 "$RADII_FILE" | while read -r znucl n rms_fm; do
    # Skip empty lines or lines that don't have the expected format
    if [ -z "$znucl" ] || [ -z "$rms_fm" ]; then
        continue
    fi

    # Clean up potential extra spaces or tabs in rms_fm (important for Fortran read)
    rms_fm_clean=$(echo "$rms_fm" | xargs)

    echo "Processing Z = $znucl, R_fm = $rms_fm_clean"

    # Run for 2s state
    echo "  Running for state 2s..."
    printf "%s\n%s\n%s\n" "$znucl" "2s" "$rms_fm_clean" | "$EXECUTABLE"

    # Check exit status (optional, but good practice)
    if [ $? -ne 0 ]; then
        echo "  Warning: Run for Z=$znucl, state=2s failed."
    else
        echo "  Completed state 2s."
    fi

    # Run for 2p state
    echo "  Running for state 2p..."
    printf "%s\n%s\n%s\n" "$znucl" "2p" "$rms_fm_clean" | "$EXECUTABLE"

    # Check exit status
    if [ $? -ne 0 ]; then
        echo "  Warning: Run for Z=$znucl, state=2p failed."
    else
        echo "  Completed state 2p."
    fi

done

echo "All runs finished."
echo "Output files (e.g., ${znucl}_polynomial_table_2s.txt) should be located in $RUN_DIR"

exit 0