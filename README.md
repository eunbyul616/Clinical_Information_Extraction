# Clinical_Information_Extraction
It is a system that extracts information from clinical data based on rules. In particular, it aims to extract important information of the surgical pathology record of electronic health records(EHRs).

## Extracting information from surgical pathology records

### Information to extract

* organ
* location
* operation name
* histologic diagnosis
* tumor size
* number of test (#1_test ~ #13_test) 
* number of positive (#1 ~ #13) 

### Installation and Requirements

```sh
$ pip install -r requirements.txt
```

### To execute code

```sh
$ python extract_information.py
```
