[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/isabela42/GRADE2">
    <!-- <img src="images/pipelineSimplels.png" alt="Logo" width=400>-->
  </a>

  <h3 align="center">GRADE2</h3>

  <p align="center">
    General RNAseq Analysis for Differential Expression version 2
  </p>
</p>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#overview">Overview</a>
    </li>
    <li>
      <a href="#pipeline-prerequisites">Pipeline prerequisites</a>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#license">License</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## Overview

This repository contains the scripts used in our in-house RNAseq analysis pipeline that runs on a PBS cluster.

Each script file foccus on one part of the analysis, from file preparation and quality assessment to differential expression.

* Script 000: Index building (001 - Kallisto; 002 RSEM; 003 STAR)
* Script 010: Quality check raw files (011 - FastQC; 012 - MultiQC)
* Script 020: Trim reads of adapters (021 - Trimmomatic)
* Script 030: Quality check raw files (031 - FastQC; 032 - MultiQC)
* Script 040: Quantify reads (041 - Kallisto)
* Script 050: Create Kallisto count tables (051 - Kallisto)
* Script 060: Alignment (061 - STAR)
* Script 070: Process alignment (071 - SAMtools; 072 - NovoSort)
* Script 080: Quantify reads (081 - RSEM)
* Script 100: Differential Expression Analysis (101 - EdgeR)

<!-- GETTING STARTED -->
## Pipeline Prerequisites

To get a local copy up and running, make sure you have each script file prerequisites instaled and up to date.

<!-- USAGE EXAMPLES -->
## Usage

Each script can be executed in a PBS cluster by using the following command line:
 
```sh
bash script-name.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"
```

<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request


<!-- ACKNOWLEDGEMENTS
## Acknowledgements

* []()
* []()
* []() -->


<!-- CONTACT -->
## Contact

Isabela Almeida - mb.isabela42@gmail.com

Project Link: [https://github.com/isabela42/GRADE2](https://github.com/isabela42/GRADE2)

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/isabela42/NB-lncRNAs.svg?style=for-the-badge
[contributors-url]: https://github.com/isabela42/GRADE2/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/isabela42/NB-lncRNAs.svg?style=for-the-badge
[forks-url]: https://github.com/isabela42/GRADE2/network/members
[stars-shield]: https://img.shields.io/github/stars/isabela42/NB-lncRNAs.svg?style=for-the-badge
[stars-url]: https://github.com/isabela42/GRADE2/stargazers
[issues-shield]: https://img.shields.io/github/issues/isabela42/NB-lncRNAs.svg?style=for-the-badge
[issues-url]: https://github.com/isabela42/GRADE2/issues
[license-shield]: https://img.shields.io/github/license/isabela42/NB-lncRNAs.svg?style=for-the-badge
[license-url]: https://github.com/isabela42/GRADE2/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://www.linkedin.com/in/isabela42/