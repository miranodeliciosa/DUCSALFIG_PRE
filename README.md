# Tiltanic Experiment

**Tiltanic** is a MATLAB-based Psychtoolbox experiment designed to study shape perception and attentional filtering using flickering bar arrays. Participants detect target shapes formed by oriented bars embedded in a dynamic background.

---

## Overview

**Goal:**  
Participants identify specific shapes (e.g., triangle, square, bar) formed by rotating subsets of bars within a grid.

**Core Features:**
- Bar grid stimulus with flickering colors (Static Bar Array, SBA)
- Foreground motion distraction (Random Dot Kinematogram, RDK)
- Target shape detection via orientation and contrast changes
- Supports randomized trials and training mode

---

## Repository Structure

```
Tiltanic/
│
├── run_Tiltanic.m                 # Main entry point
├── rand_Tiltanic.m               # Trial and event randomization
├── pres_Tiltanic.m               # Trial presentation & response collection
├── generateBarTextures.m         # Stimulus texture creation (SBA)
├── extractColorIndices.m         # Extracts bar locations from shape images
├── generate_event_onset_continuous.m # Background/foreground event jittering
├── /images/                      # Shape mask images (.tiff)
└── README.md                     # This file
```

---

## Getting Started

### Prerequisites

- MATLAB R2020b or later
- [Psychtoolbox 3](http://psychtoolbox.org/)
- Working directory with all `.m` files and `/images/` folder

### Run the Experiment

```matlab
run_Tiltanic
```

You’ll be prompted to enter experiment parameters such as participant ID, run mode (training/experiment), etc.

---

## Shape Image Setup

All target shapes are defined as `.tiff` files in `/images/`.

**Requirements:**
- Format: TIFF
- Dimensions match SBA grid size (e.g., 12×12)
- Dark pixels (RGB ≤ 180) define the shape
- Each image represents one shape class

---

## Function Summary

| Function                      | Description |
|------------------------------|-------------|
| `run_Tiltanic`               | Initializes experiment, handles blocks & saving |
| `rand_Tiltanic`              | Randomizes trial conditions & event timing |
| `pres_Tiltanic`              | Presents stimuli and collects participant responses |
| `generateBarTextures`        | Generates SBA textures, positions, and orientation sequences |
| `extractColorIndices`        | Reads TIFF mask images to identify shape bar indices |
| `generate_event_onset_continuous` | Jitters event onset timings with SOA constraints |

---

## Architecture Diagram

```
run_Tiltanic
   ├── rand_Tiltanic
   │     └── generate_event_onset_continuous
   └── pres_Tiltanic
         └── generateBarTextures
               └── extractColorIndices
```

---

## Training Mode

To run the experiment in training mode, set the `flag_training` input to `1` when calling `rand_Tiltanic` or `pres_Tiltanic`. This activates different contrast levels, event counts, and trial durations.

---

## Output & Response Data

Participant responses and trial timings are stored in the `resp` and `timing` structures:
- Button press timestamps
- Correct/incorrect classification
- Reaction times
- Event metadata (contrast, shape, type)

---

## License

MIT License (see `LICENSE` file if provided)

---

## Author

**Sebastian Wehle**  
University of Leipzig, 2025  
[GitHub Profile](https://github.com/miranodeliciosa) 

---