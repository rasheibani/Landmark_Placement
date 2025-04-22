# Landmark Placement Optimization in Floorplans

This project focuses on the strategic placement of landmarks in floorplans, utilizing both real-world and synthetic data. The goal is to develop computational methods for optimal landmark allocation to aid in wayfinding and spatial navigation.

## Project Structure

- **main.py**: The main script that coordinates the landmark placement process.
- **para.py**: Contains parallel processing configurations used throughout the project.
- **computational_geometry_functions.py**: Implements various geometric functions essential for landmark computation and placement.
- **Complexity.xlsx**: Spreadsheet with data on complexity measures for the floorplans.
- **data_realworld/**: Contains real-world floorplan data.
- **data_synthetic/**: Contains synthetic floorplan data generated for testing and validation.
- **Geojson/**: Directory for storing GeoJSON files related to the floorplans.
- **results/**: Directory for saving results, such as the location and analysis of placed landmarks.
- **sub_functions/**: Contains helper functions used across different scripts.

## Requirements

To run this project, you'll need:

- Python 3.x
- Dependencies (you may list the specific libraries required, if any)

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/landmark_placement.git
   cd landmark_placement
   ```


## Usage

1. Set the parameters in `para.py` according to your project requirements.
2. Run the main script:
   ```bash
   python main.py
   ```

This will execute the landmark placement algorithms using the configurations specified in `para.py`.

## Results

The results, including the coordinates and analysis of landmark placements, are saved in the `results/` directory. These can be visualized or analyzed further as needed.

## Data

- `data_realworld/`: Contains real-world floorplans.
- `data_synthetic/`: Contains synthetic floorplans generated for testing purposes.
- `Geojson/`: Stores GeoJSON files, which are used for spatial analysis and visualization.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## References

If you use this project in your research, please consider citing the following paper:

- Arabsheibani, R., Haunert, J.-H., Winter, S., & Tomko, M. (2024). *Strategic allocation of landmarks to reduce uncertainty in indoor navigation*. Computers, Environment and Urban Systems, 114, 102198. [https://doi.org/10.1016/j.compenvurbsys.2024.102198](https://doi.org/10.1016/j.compenvurbsys.2024.102198)

  **BibTeX:**
  ```bibtex
  @article{ARABSHEIBANI2024102198,
    title = {Strategic allocation of landmarks to reduce uncertainty in indoor navigation},
    journal = {Computers, Environment and Urban Systems},
    volume = {114},
    pages = {102198},
    year = {2024},
    issn = {0198-9715},
    doi = {https://doi.org/10.1016/j.compenvurbsys.2024.102198},
    url = {https://www.sciencedirect.com/science/article/pii/S0198971524001273},
    author = {Reza Arabsheibani and Jan-Henrik Haunert and Stephan Winter and Martin Tomko},
    keywords = {Indoor navigation, Allocation, Landmark, Decision points, Route instructions, Optimization},
    abstract = {Indoor navigation systems often rely on verbal, turn-based route instructions. These can, at times, be ambiguous at complex decision points with multiple paths intersecting under angles that are not well distinguished by the turn grammar used. Landmarks can be included into turn instructions to reduce this ambiguity. Here, we propose an approach to optimize landmark allocation to improve the clarity of route instructions. This study assumes that landmark locations are constrained to a pre-determined set of slots. We select a minimum-size subset of the set of all slots and allocate it with landmarks, such that the navigation ambiguity is resolved. Our methodology leverages computational geometric analysis, graph algorithms, and optimization formulations to strategically incorporate landmarks into indoor route instructions. We propose a method to optimize landmark allocation in indoor navigation guidance systems, improving the clarity of route instructions at complex decision points that are inadequately served by turn-based instructions alone.}
  }
  ```
