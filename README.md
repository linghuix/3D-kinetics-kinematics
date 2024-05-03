# Biomechanical Analysis MATLAB Code

## Overview

This MATLAB script is designed to perform biomechanical analysis related to joint angles, angular velocities, accelerations, ground forces, and moments. The script is structured into several sections, each serving a specific purpose within the biomechanical analysis pipeline.

## Usage

1. **Function: test()**
   - This function serves as the entry point to the biomechanical analysis.
   - It initializes the workspace, sets up parameters such as time step, and calls various sub-functions to perform specific calculations.
2. **Joint Angle Calculation**
   - The `get_joint_angle()` function calculates joint angles based on the input data.
   - The `plot_joint_angle()` function visualizes the joint angles over time.
3. **Angular Velocity and Acceleration Calculation for Joints**
   - Angular velocities and accelerations for different joints are calculated using numerical differentiation.
4. **Acceleration Calculation**
   - Bone coordinates are loaded, and absolute acceleration for each bone is computed.
5. **Ground Force Calculation**
   - Ground forces and their positions are determined.
6. **Moment Calculation**
   - Moments are calculated based on forces, moments of inertia, and other biomechanical parameters.

## Dependencies

- MATLAB R2022 or later

## How to Run

1. Ensure MATLAB is installed on your system.
2. Open MATLAB and navigate to the directory containing the script.
3. Run the script by typing `test()` in the MATLAB command window.

## Contributing

Contributions to this project are welcome! If you have suggestions for improvements or new features, feel free to open an issue or submit a pull request.

## License

This project is licensed under the MIT License.
