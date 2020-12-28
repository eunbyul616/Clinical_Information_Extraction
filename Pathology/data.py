import pandas as pd
import os

class Dataset:
    """

    load data

    """
    def __init__(self, path, file_name):
        self.path = path
        self.file_name = file_name
    
    def load_data(self, categories, columns_=None):
        self.df = pd.DataFrame()
        
        for category in categories:
            tmp = pd.read_excel(os.path.join(self.path, self.file_name), sheet_name = category)
            if columns_ != None:
                tmp = tmp[columns_]
            
            tmp['category'] = category
            print('{} shape: {}'.format(category, tmp.shape))
            self.df = pd.concat([self.df, tmp])
    
    def __str__(self):
        return 'total shape: {}'.format(self.df.shape)