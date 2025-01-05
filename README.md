# Photon Trajectory Simulation Project

This project simulates photon trajectories around a Schwarzschild black hole. It calculates the path of photons and visualizes their journey toward the event horizon and the accretion disk.

## Table of Contents

- [Project Overview](#project-overview)
- [Installation](#installation)
- [Usage](#usage)
- [File Descriptions](#file-descriptions)
- [Contributing](#contributing)

## Project Overview

This simulation allows you to study the behavior of photons as they travel around a Schwarzschild black hole, showcasing the bending of light in the strong gravitational field. The main objective is to trace multiple photon paths originating from different angles and observe how they interact with the black hole's event horizon and accretion disk.

The project utilizes numerical methods to solve the equations of motion for photons and provides visualizations of photon trajectories. These visualizations can be analyzed to understand phenomena like gravitational lensing and photon ring formation.

## Installation

To run this project, you'll need Python and a few libraries. Follow these steps to set everything up:

### Prerequisites

- Ensure you have Python 3.x installed.
- You'll need the following Python libraries: NumPy and Matplotlib.

### Installation Steps

1. Clone this repository to your local machine.

2. (Optional but recommended) Create a virtual environment to isolate the project dependencies.

3. Install the required dependencies. You can do this by installing the necessary libraries through your package manager, or if using a virtual environment, activate it and install the libraries in the environment.

4. If you are using a `requirements.txt` file, simply install the dependencies listed in the file.

## Usage

1. Ensure all dependencies are installed as described in the [Installation](#installation) section.
   
2. To run the photon trajectory simulation, execute the main script:

   - `MultipleRayPaths.py` — This script generates photon paths from different initial angles and visualizes them around the black hole.

3. The program will output a visualization of photon paths, along with information about the photon’s interaction with the accretion disk.

4. You can modify the initial parameters in the scripts to simulate different scenarios, such as changing the radius from which the photons are emitted or adjusting the number of rays.

## File Descriptions

- **MultipleRayPaths.py**: This is the main script for simulating multiple photon trajectories around a Schwarzschild black hole. It generates photon paths from different initial angles and visualizes them in 2D space.

- **SingleRayPath.py**: Contains the function to compute the path of a single photon. It is used in the main script to compute individual trajectories.

- **coloredLine.py**: Provides the functionality to plot photon paths with color coding based on the time it takes for photons to reach the accretion disk.

- **cubicequationsolver.py**: Implements a cubic equation solver used in solving certain equations related to photon trajectories.

## Contributing

Contributions to this project are welcome. If you'd like to contribute, please fork the repository, make your changes, and submit a pull request with a clear explanation of your modifications. If you encounter any issues or have suggestions, feel free to open an issue.


