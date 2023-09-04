# Gene Planet Assignment

## Description

The Gene Planet Assignment project is designed to analyze genomic data from BAM files, calculate average coverage across the genome, and generate reports summarizing the analysis results.

## How to Build Docker Image

To build the Docker image for this project, follow these steps:

1. Clone the repository to your local machine:
git clone https://github.com/AljazMekinda/gene-planet-homework.git

2. Navigate to the project directory:

3. Build the Docker image using the provided Dockerfile: docker build -t gene-planet:0.3 .

## How to Run and Mount Volumes

Once you have built the Docker image, you can run the project and mount volumes to provide input data and retrieve output reports. To run the project with volume mounting:

```bash
docker run -v /path/to/data:/data -v /path/to/reports:/reports -v /path/to/runner_config:/config/runner_config gene-planet:0.3  /bin/bash

```

Make sure to replace **/path/to/data** and **/path/to/reports** with the local directory containing your BAM files and **/path/to/runner_config** with the directory containing your configuration file.

## Config file

```json
{
  "project_name": "gene-planet-assignment",
  "bam_file": "data/79734fd3-87c2-4a79-a25d-4c4465c432e2_alignment-bam.bam",
  "report_path": "reports/genome_report.pdf"
}
```

You can customize the project_name, bam_file, and report_path as needed to specify your project details and file locations.


## Requirements

    Docker
