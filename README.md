#### ABacAPP: Anti-BACterial Agent Potency Predictor

![Logo](https://github.com/naeemmrz/ABacAPP/blob/main/logo.png?raw=true)

<div style="text-align: justify"> 
ABacAPP is a machine learning-based web application capable of predicting the antibacterial potency of any compound against β-lactamases (inlcluding β-lactamase, AmpC, Bla2, KPC-2, TEM-1, BRO-1, Class C, Class D, Type II, L1, NDM-1, SHV-1, SHV-5, OXA-1, and Metallo β-lactamase 1 and 2), DNA Gyrases (inlcluding DNA Gyrase, GyraseA, and GyraseB), Bacterial Dihydrofolate Reductases, and Penicillin-Binding Proteins (inlcuding 1B, 2, 2B, 2X, MecA, and D-alanyl-D-alanine carboxypeptidase).
</div>


## Authors
**Naeem Abdul Ghafoor¹²³** & **Omur Baysal³**
###### ¹Universidad Autónoma de Barcelona, 08193 Cerdanyola, Spain.
###### ²Faculty of Biology, University of Barcelona, 08028 Barcelona, Spain
###### ³Mugla Sitki Kocman University, Faculty of Science, Department of Molecular Biology and Genetics, Mugla, Turkey.
###### [Click Here For The Corresponding Article](https://)


## Usage
## [![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://naeemmrz-abacapp-abacapp-qpjwr6.streamlitapp.com)
- Click on the "Open in Streamlit" badge above. 
- Enter the SMILE for the compound of your interest.
- Predicted IC50 will be printed out within a few seconds.
  
  
## Run Locally
Download the project

```bash
  wget https://github.com/naeemmrz/ABacAPP.git
```

Unzip and enter the project directory

```bash
  unzip ABacAPP-main.zip
  cd ABacAPP-main
```

Install dependencies

```bash
  pip install -r requirements.txt  #System level packages are provided in packages.txt
```

Start the application

```bash
  streamlit run MDM2pred.py
```

The Application will open in your default browser with the same interface as the online version.


## Acknowledgements
<div style="text-align: justify"> 
ABacAPP is developed and maintained by Naeem A. from <a href="http://ynlab.mu.edu.tr/">YNLab</a> at the Department of Molecular Biology and Genetics, Mugla Sitki Kocman University, Mugla, Turkey. The development of ABacAPP was funded by the The Scientific and Technological Research Institution of Turkey (TUBITAK) under the project number 122E082 and supervised by Prof. Dr. Omur Baysal from the Molecular Microbiology Laboratory at the Department of Molecular Biology and Genetics, Faculty of Science, Mugla Sitki Kocman University, Mugla, Turkey. The numerical calculations that lead to the development of ABacAPP were partially performed at TUBITAK ULAKBIM, High Performance and Grid Computing Center (TRUBA resources). The authors/developers are grateful to TUBITAK and TRUBA for the funding and resources they've provided.
</div>
