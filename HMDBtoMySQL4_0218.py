# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 13:01:06 2018

@author: LiWenHuang
"""
#連接mysql、建立表格與欄位
import mysql.connector as mydb 
conn = mydb.connect(host = '127.0.0.1',user = 'root',passwd = '1234',
    db = 'compound',port = 3306,charset = 'utf8')
cursor = conn.cursor()
tables = '''create table hmdbdata20190215(creation_date text not null,update_date text not null,
            accession text not null,name text not null,description text not null,
            synonyms_synonym text not null,chemical_formula text not null,average_molecular_weight text not null,
            monisotopic_molecular_weight text not null,iupac_name text not null,
            traditional_iupac varchar(1000) not null,cas_registry_number text not null,
            smiles text not null, inchi text not null, inchikey text not null,
            status text not null,ontology_origin text ,
            ontology_biofunction text ,cellular_location text not null,
            spectra_type text not null,spectra_spectrum_id text not null,
            biospecimen_locations text not null,tissue_locations_tissue text not null,
            pathways_name longtext not null,pathways_smpdb_id longtext not null,
            pathways_kegg_map_id longtext not null,disease_name text not null,
            disease_omim_id text not null,drugbank_id text not null,foodb_id text not null,
            knapsack_id text not null,chemspider_id text not null,kegg_id text not null,
            biocyc_id text not null,bigg_id text not null,wikipidia text not null,
            nugowiki text not null,metagene text not null,metlin_id text not null,
            pubchem_compound_id text not null,het_id text not null,chebi_id text not null)'''
cursor.execute(tables)
print("table finish")
#取出標籤內容存成一個list
import xml.etree.ElementTree as ET
tree = ET.parse('hmdb_metabolites.xml')
root = tree.getroot()
label = '{http://www.hmdb.ca}'
label1 = ['creation_date','update_date','accession','name','description','chemical_formula','average_molecular_weight',
         'monisotopic_molecular_weight','iupac_name','traditional_iupac','cas_registry_number','smiles','inchi','inchikey',
         'status','drugbank_id','foodb_id','knapsack_id','chemspider_id','kegg_id','biocyc_id','bigg_id','wikipidia','nugowiki',
         'metagene','metlin_id','pubchem_compound_id','het_id','chebi_id']
label2 = ['synonyms','cellular_locations','biospecimen_locations','tissue_locations']
label2_1 =['synonym','cellular','biospecimen','tissue']
label3 = ['spectra','spectra','pathways','pathways','pathways','diseases','diseases']
label3_1 = ['spectrum','spectrum','pathway','pathway','pathway','disease','disease']
label3_1_1 = ['type','spectrum_id','name','smpdb_id','kegg_map_id','name','omim_id']
for a in range(len(label1)):
    label1[a] = label + label1[a] 
for b in range(len(label2)):
    label2[b] = label + label2[b]
    label2_1[b] = label + label2_1[b]
for c in range(len(label3)):
    label3[c] = label + label3[c] 
    label3_1[c] = label + label3_1[c] 
for d in range(len(label3_1_1)):
    label3_1_1[d] = label + label3_1_1[d] 
data1 = []
data2 = []
data3 = []
data4 = []
count = 0
for subroot in root:
    count +=1    
    for i in range(0,len(label1),1):
        if len(subroot.findall(label1[i])) == 0 and subroot.text.strip()=='' :
            text="no exist"
            data1.append(text)
        else:
            for child in subroot.findall(label1[i]):
                if len(child)==0 and child.text==None:
                    text='None'
                    data1.append(text)
                elif len(child) == 0 and child.text.strip()=='':
                    text="NULL"
                    data1.append(text)
                elif len(child)==0 and child.text.strip():
                    text=child.text
                    data1.append(text)
    for j in range(0,len(label2),1):
        if len(subroot.findall(label2[j])) == 0 and subroot.text.strip()=='' :
            text="no exist"
            data2.append(text)
        else:
            for child in subroot.findall(label2[j]):
                if len(child)==0 and child.text==None:
                    text='None'
                    data2.append(text)
                elif len(child) == 0 and child.text.strip()=='':
                    text="NULL"
                    data2.append(text)
                elif len(child)==0 and child.text.strip():
                    text=child.text
                    data2.append(text)
                else:
                    subdata=[]
                    for subchild in child.findall(label2_1[j]):
                        text=subchild.text
                        subdata.append(text)
                    combine_subdata = ";".join(subdata)
                    data2.append(combine_subdata)
    for k in range(0,len(label3),1):
        data3.clear()
        if len(subroot.findall(label3[k])) == 0 and subroot.text.strip()=='' :
            text="no exist"
            data3.append(text)
        else:
            for child in subroot.findall(label3[k]):
                if len(child)==0 and child.text==None:
                    text='None'
                    data3.append(text)
                elif len(child) == 0 and child.text.strip()=='':
                    text="NULL"
                    data3.append(text)
                elif len(child)==0 and child.text.strip():
                    text=child.text
                    data3.append(text)
                else:
                    for subchild in child.findall(label3_1[k]):
                        subsubdata=[]
                        for subsubchild in subchild.findall(label3_1_1[k]):
                            if subsubchild.text==None:
                                text="None"
                                subsubdata.append(text)
                            elif subsubchild.text.strip()=='':
                                text="NULL"
                                subsubdata.append(text)
                            else:
                                text = subsubchild.text
                                subsubdata.append(text)
                        combine_subdata2 = ";".join(subsubdata)
                        data3.append(combine_subdata2)
        combine_data = ';'.join(data3)
        data4.append(combine_data)
    data1.insert(5,data2[0])
    data1.insert(16,'')
    data1.insert(17,'')
    data1.insert(18,data2[1])
    data1.insert(19,data4[0])
    data1.insert(20,data4[1])
    data1.insert(21,data2[2])
    data1.insert(22,data2[3])   
    data1.insert(23,data4[2])
    data1.insert(24,data4[3])
    data1.insert(25,data4[4])
    data1.insert(26,data4[5])
    data1.insert(27,data4[6])
#將資料所有的值去逗號
    for t in range(len(data1)):
        data1[t] = data1[t].replace(',','*')
        data1[t] = data1[t].replace('α','alpha')
        data1[t] = data1[t].replace('β','beta')
        data1[t] = data1[t].replace('ω','omega')
#將list中資料一筆一筆丟入mysql 
    cursor.execute("insert into hmdbdata20190215 values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\
                                            %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",[str(data1[p]) for p in range(len(data1))])
    data1.clear()
    data2.clear()
    data3.clear()
    data4.clear()
    subdata.clear()
    subsubdata.clear()
    print(count)
cursor.close()
print('datafinish') 