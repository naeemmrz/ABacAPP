import warnings
warnings.filterwarnings("ignore")

import os
import PIL
import pickle
import pandas as pd
import streamlit as st
import streamlit_option_menu
from streamlit_option_menu import option_menu

st.set_page_config(page_title='ABacAPP', page_icon=PIL.Image.open('logo.png'))

selected = option_menu(menu_title=None,
			options=["Home", "Run ABacAPP", "Batch Mode", "About"],
			icons=["house", "book", "clock", "envelope"],
			menu_icon="cast",
			default_index=0,
			orientation="horizontal")
			
if selected == "Home":
	col1, col2, col3 = st.columns([0.1,8,0.1])
	with col1:
		st.write("")
	with col2:
		image = PIL.Image.open('logo.png')
		st.image(image, caption=None, use_column_width=True, width=500)
	with col3:
		st.write("")
	st.markdown("<h4 style='text-align: center; color: white;'> ABacAPP: Anti-BACterial Agent Potency Predictor</h4>", unsafe_allow_html=True)
	st.write("\n")
	st.write("\n")
	
	st.markdown("<h6 style='text-align: justify; color: #C0C0C0;'> ABacAPP is a machine learning-based web application capable of predicting the antibacterial potency of any compound against β-lactamases (including β-lactamase, AmpC, Bla2, KPC-2, TEM-1, BRO-1, Class C, Class D, Type II, L1, NDM-1, SHV-1, SHV-5, OXA-1, and Metallo β-lactamase 1 and 2), DNA Gyrases (including DNA Gyrase, GyraseA, and GyraseB), Bacterial Dihydrofolate Reductases, and Penicillin-Binding Proteins (including 1B, 2, 2B, 2X, MecA, and D-alanyl-D-alanine carboxypeptidase). </h6>", unsafe_allow_html=True)
	st.write("\n")
	
	st.markdown("<h6 style='text-align: justify; color: #C0C0C0;'> To use the app, click the ''Run ABacAPP'' tab on the top of the page and enter the SMILES representation of the target compound in the appropriate field. Futhermore, to run prediction for several compounds at one, click the ''Batch Mode'' tab and upload CSV file with smiles in the first column and identifiers (arbitrary but necessary to run the application) in the second column (without headers or index column), a template is available for download in the respective tab.\n </h6>", unsafe_allow_html=True)
	st.write("\n")
	
	st.markdown("<h6 style='text-align: justify; color: #C0C0C0;'> The β-lactamases model is based on the Random Forest regressor and the Klekota Roth Count features (shape=4860), the model was trained on 2045 datapoints and achives R², MSE, RMSE, and MAE of 0.87, 0.42, 0.65, and 0.44 respectively.\n The DNA Gyrases model is based on the K-Nearest Neighbors regressor and the Circular Fingerprint features (radius=2 and shape=1024), the model was trained on 813 datapoints and achives R², MSE, RMSE, and MAE of 0.82, 0.67, 0.80, and 0.36 respectively.\n The Bacterial Dihydrofolate Reductase model is based on the Random Forest regressor and the Klekota Roth Count features (shape=4860), the model was trained on 306 datapoints and achives R², MSE, RMSE, and MAE of 0.81, 0.74, 0.82, and 0.54 respectively.\n The Penicillin-Binding Proteins model is based on the Extra Trees regressor and the Substructure Fingerprint Count features (shape=307), the model was trained on 132 datapoints and achives R², MSE, RMSE, and MAE of 0.88, 0.22, 0.44, and 0.24 respectively (All benchmarks are for LOOCV on the test set with CV=10).\n </h6>", unsafe_allow_html=True)
	st.write("\n")
	st.write("\n")
	st.write("\n")
	st.markdown("<h6 style='text-align: justify; color: white;'> Refer to **article link** for more details. Please cite ABacAPP as ''''pending'''' if you have used it in your research and/or provide credits to the orginal App/authors if used for any other purpose.\n </h6>", unsafe_allow_html=True)

	
if selected == "Run ABacAPP":
	st.title(f"Welcome to ABacAPP")
	user_input = st.text_input("", "CC(=O)NC1=CC=C(C=C1)O")
	if user_input is None:
		st.write(f"Waiting for user input")
	else:
		smile = user_input
	try:
		from rdkit import Chem
		m = Chem.MolFromSmiles(smile)
		csmi = Chem.rdmolfiles.MolToSmiles(m)
	except:
		st.write(f"Please provide a valid SMILE")
		st.stop()
		
	def single_KrF(csmi):
		df = pd.DataFrame({'smiles': csmi, 'id': 'target'}, index=[0])
		df.to_csv('molecule.smi', sep='\t', header=False, index=False)
		os.system('java -Xms1G -Xmx1G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/KlekotaRothFingerprintCount.xml -dir ./molecule.smi -file descriptors_output.csv')
		df = pd.read_csv('descriptors_output.csv', sep=',')
		os.system('rm -rf descriptors_output.csv molecule.smi')
		df = df.drop('Name', axis=1)
		return df

	def single_SbF(csmi):
		df = pd.DataFrame({'smiles': csmi, 'id': 'target'}, index=[0])
		df.to_csv('molecule.smi', sep='\t', header=False, index=False)
		os.system('java -Xms1G -Xmx1G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/SubstructureFingerprintCount.xml -dir ./molecule.smi -file descriptors_output.csv')
		df = pd.read_csv('descriptors_output.csv', sep=',')
		os.system('rm -rf descriptors_output.csv molecule.smi')
		df = df.drop('Name', axis=1)
		return df

	def single_ExC2R1024(csmi):
		import deepchem
		featurizer = deepchem.feat.CircularFingerprint(size=1024, radius=2, chiral=True)
		features = featurizer.featurize(csmi)
		return features
	def pIC50_2_mM(pIC50):
		from math import e
		ic50_mM = float( e ** (-pIC50)) * (10 ** 6) / 1000
		return ic50_mM
		
	BL_Model = pickle.load(open('models/BL_KrF_RF.pkl', 'rb'))
	DR_Model = pickle.load(open('models/DR_KrF_RF.pkl', 'rb'))
	GY_Model = pickle.load(open('models/GY_ExC2R1024_KNN.pkl', 'rb'))
	NPN_Model = pickle.load(open('models/NPN_SbF_ETR.pkl', 'rb'))
	try:
		BL_pIC50 = BL_Model.predict(single_KrF(csmi))
		DR_pIC50 = DR_Model.predict(single_KrF(csmi))
		GY_pIC50 = GY_Model.predict(single_ExC2R1024(csmi))
		NPN_pIC50 = NPN_Model.predict(single_SbF(csmi))
	except:
		st.write(f"Oops!, counldn't calculate features for the provided SMILE, please recheck it or try a different one")
		st.stop()
	st.write('\n')
	st.write('\n')	
	st.write(f"Predicted IC50 against β-lactamase:                      {round(pIC50_2_mM(BL_pIC50), 3)} mM")
	st.write(f"Predicted IC50 against acterial Dihydrofolate Reductase: {round(pIC50_2_mM(DR_pIC50), 3)} mM")
	st.write(f"Predicted IC50 against DNA Gyrase:                       {round(pIC50_2_mM(GY_pIC50), 3)} mM")
	st.write(f"Predicted IC50 against Penicillin-Binding Proteine:      {round(pIC50_2_mM(NPN_pIC50), 3)} mM")
	
if selected == "Batch Mode":
	st.title(f"Welcome to ABacAPP - Batch Mode")
	uploaded_file = st.file_uploader("Upload your input CSV file", type=["csv"])
	def smiles2canonical(df):
		from rdkit import Chem
		canon_smiles = []
		smiles = df['smiles']
		for i in range(len(smiles)):
			m = Chem.MolFromSmiles(smiles[i])
			try:
				csmi = Chem.rdmolfiles.MolToSmiles(m)
			except:
				csmi = 'uncanonicalizable'
			canon_smiles.append(csmi)
		df['canonical_smiles'] = canon_smiles
		df = df.loc[df['canonical_smiles'] != 'uncanonicalizable']
		df.reset_index(inplace=True, drop=True)
		df = df.drop('smiles', axis=1)
		return df 

	def batch_KrF(df):
		df[['canonical_smiles', 'ids']].to_csv('molecule.smi', sep='\t', header=False, index=False)
		os.system('java -Xms1G -Xmx1G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/KlekotaRothFingerprintCount.xml -dir ./molecule.smi -file descriptors_output.csv')
		KrF = pd.read_csv('descriptors_output.csv', sep=',')
		KrF = KrF.drop('Name', axis=1)
		os.system('rm -rf descriptors_output.csv molecule.smi')
		return KrF

	def batch_SbF(df):
		df[['canonical_smiles', 'ids']].to_csv('molecule.smi', sep='\t', header=False, index=False)
		os.system('java -Xms1G -Xmx1G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/SubstructureFingerprintCount.xml -dir ./molecule.smi -file descriptors_output.csv')
		SbF = pd.read_csv('descriptors_output.csv', sep=',')
		SbF = SbF.drop('Name', axis=1)
		os.system('rm -rf descriptors_output.csv molecule.smi')
		return SbF

	def batch_ExC2R1024(df):
		import deepchem
		smiles = df['canonical_smiles']
		featurizer = deepchem.feat.CircularFingerprint(size=1024, radius=2, chiral=True)
		ExC2R1024 = featurizer.featurize(smiles)
		return pd.DataFrame(ExC2R1024)

	def pIC50_2_mM(pIC50):
		from math import e
		ic50_mM = float( e ** (-pIC50)) * (10 ** 6) / 1000
		return ic50_mM

	def batch_predictions(df):
		df = smiles2canonical(df)
		st.write('Featurizing ...')
		KrF = batch_KrF(df)
		SbF = batch_SbF(df)
		ExC2R1024 = batch_ExC2R1024(df)
		st.write('Loading Models ..')
		BL_Model = pickle.load(open('models/BL_KrF_RF.pkl', 'rb'))
		DR_Model = pickle.load(open('models/DR_KrF_RF.pkl', 'rb'))
		GY_Model = pickle.load(open('models/GY_ExC2R1024_KNN.pkl', 'rb'))
		NPN_Model = pickle.load(open('models/NPN_SbF_ETR.pkl', 'rb'))
		st.write('Making Predictions ...')
		BL_pIC50 = BL_Model.predict(KrF)
		DR_pIC50 = DR_Model.predict(KrF)
		GY_pIC50 = GY_Model.predict(ExC2R1024)
		NPN_pIC50 = NPN_Model.predict(SbF)
		st.write('Preparing Results ...')
		df['BL_pIC50'] = BL_pIC50
		df['DR_pIC50'] = DR_pIC50  
		df['GY_pIC50'] = GY_pIC50
		df['NPN_pIC50'] = NPN_pIC50

		df['Beta-Lactamase IC50 mM'] = df['BL_pIC50'].apply(pIC50_2_mM).round(decimals=3)
		df['Dihydrofolate Reductase IC50 mM'] = df['DR_pIC50'].apply(pIC50_2_mM).round(decimals=3)
		df['DNA Gyrase IC50 mM'] = df['GY_pIC50'].apply(pIC50_2_mM).round(decimals=3)
		df['Penicillin-Binding Protein IC50 mM'] = df['NPN_pIC50'].apply(pIC50_2_mM).round(decimals=3)
	  
		DF = df.drop(['BL_pIC50', 'DR_pIC50', 'GY_pIC50', 'NPN_pIC50'], axis=1)
		return DF

	@st.cache
	def convert_df(df):
		return df.to_csv(index=False, header=None).encode('utf-8')
			
	if uploaded_file is None:
		sample = pd.read_csv('sample_batch_input.csv', sep=',', header=None)
		sample_csv = convert_df(sample)
		col1, col2, col3 = st.columns([2.5,8,1])
		with col1:
			pass
		with col2:
			st.download_button('Click here to download a template for the batch mode', sample_csv, mime='text/csv', file_name='sample_batch_input.csv', key='download-csv')
		with col3:
			pass
		st.write(f"Waiting for user input")
	else:
		st.write(f"Processing user input ... ")
		input_csv = uploaded_file
		df = pd.read_csv(input_csv, sep=',', header=None, names=['smiles', 'ids'])
		try:
			DF = batch_predictions(df)
			st.write(DF)
			csv = convert_df(DF)
		except:
			st.write(f"Oops, counld processes the input, please follow the given template for batch mode.")
			st.stop()
		st.download_button('Download results as csv', csv, mime='text/csv', file_name='ABacAPP_batch_predictions.csv', key='download-csv')

if selected == "About":
	st.write("\n")
	st.write("\n")
	st.markdown("""<h6 style='display: block; text-align: justify; color: #C0C0C0;'> ABacAPP is developed and maintained by Naeem A. from <a href="http://ynlab.mu.edu.tr/">YNLab</a> at the Department of Molecular Biology and Genetics, Mugla Sitki Kocman University, Mugla, Turkey. </h6>""", unsafe_allow_html=True)
	st.write('\n')
	
	st.markdown("""<h6 style='text-align: justify; color: #C0C0C0;'> The development of ABacAPP was funded by the The Scientific and Technological Research Institution of Turkey (TUBITAK) under the project number 122E082 and supervised by <a href="https://www.mu.edu.tr/tr/personel/omurbaysal">Prof. Dr. Omur Baysal</a> from the Molecular Microbiology Laboratory at the Department of Molecular Biology and Genetics, Faculty of Science, Mugla Sitki Kocman University, Mugla, Turkey. The numerical calculations that lead to the development of ABacAPP were partially performed at TUBITAK ULAKBIM, High Performance and Grid Computing Center (TRUBA resources). The authors/developers are grateful to TUBITAK and TRUBA for the funding and resources they've provided.</h6>""", unsafe_allow_html=True)
	st.write('\n')

	st.markdown("<h6 style='text-align: justify; color: #C0C0C0;'> ABacAPP is licensed under the MIT License Copyright (c) 2022 Naeem A., Email to merzanaeem007 [at] gmail.com to report any bugs or recommend further features/improvments. </h6>", unsafe_allow_html=True)
	st.write('\n')
	st.write('\n')
	
	st.markdown("<h6 style='text-align: justify; color: white;'> Refer to **article link** for more details. Please cite ABacAPP as ''''pending'''' if you have used it in your research and/or provide credits to the orginal App/authors if used for any other purpose.\n </h6>", unsafe_allow_html=True)
