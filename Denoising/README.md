# Image Denoising Project

## Project Overview
This repository is dedicated to exploring LP relaxations for image denoising, particularly focusing on synthetic matrices and QR codes. The project is structured into two main components: CRAMA and QR_train_set, each containing specific types of data used for denoising experiments.

## Repository Structure
- **CRAMA/**: Contains synthetic matrices for three types of images:
  - **CEN**: Center rectangle.
  - **TL**: Top-left rectangle.
  - **CROSS**:

  Each pattern type is available in two sizes: 15x15 and 100x100 pixels.

- **QR_train_set/**: A collection of 50x50 QR codes, each encoded with random information. These QR codes are used to learn potential values for image denoising tasks.

- **qr_200_100.txt**: A 200x200 QR code used for testing the denoising methods.

- **CRAMA_denoising.ipynb**: An Jupyter Notebook that includes all necessary steps to reproduce the experiments on synthetic images.

- **QR_code_denoising.ipynb**: An Jupyter Notebook that includes all necessary steps to reproduce the experiments on QR code denoising.