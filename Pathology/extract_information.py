import pandas as pd
import numpy as np
import os
import glob
import re
from tqdm import tqdm
import collections
import ast
import json
from util import *
from data import Dataset
tqdm.pandas()

"""

* set variables
    DATA_PATH: 데이터 경로
    RESULT_PATH: 최종 결과 파일 저장 경로
    result_fname: 최종 결과 파일명
    file_name: 데이터 파일명
    categories: xlsx 파일 sheet명
    columns_: 데이터 파일 내 column명
    
* 최종 결과 파일
    - final_result.csv 
    - result_fname
    - null.csv : 추출 불가 행
    
"""

DATA_PATH = './data/'
RESULT_PATH = './result/'
os.makedirs(RESULT_PATH, exist_ok=True)
result_fname = 'final_result_2020.csv'

# file_name = 'selected_outcome1.xlsx'
# categories = ['other_needle', 'lung_excision', 'lung_needle', 'washing_fluid', 'lymph_node1', 'lymph_node2']
# columns_ = ['IRB', 'idx', 'code', 'code_NM', 'day', 'main_script']

file_name = '수술병리_2020.xlsx'
categories = ['Sheet1']
columns_ = None

if __name__ == '__main__':
    # load data
    print('START LOAD DATA ...')
    dataset = Dataset(DATA_PATH, file_name)
    dataset.load_data(categories, columns_)
    print(dataset)

    dataset.df.columns = ['IRB', 'idx', 'day', 'code', 'code_NM', 'main_script', 'category']
    print(dataset.df.columns)
    print(dataset.df.head())

    columns_ = list(dataset.df.columns.values)

    print('START EXTRACT ...')
    tmp_df = dataset.df.copy()
    tmp_df.dropna(inplace=True)
    tmp_df.reset_index(inplace=True, drop=True)
    tmp_df['extract_diagnosis'] = tmp_df['main_script'].progress_apply(lambda x: extract_diagnosis(x))
    tmp_df['split_sent'] = tmp_df['extract_diagnosis'].progress_apply(lambda x: split_sent(x))

    rows = []
    for idx in tqdm(range(tmp_df.shape[0])):
        for info in tmp_df['split_sent'].iloc[idx]:
            info, histology = info[0], info[1]
            list_ = tmp_df[columns_ + ['extract_diagnosis']].iloc[idx].tolist()
            list_.extend([info, histology])
            rows.append(list_)

    result_df = pd.DataFrame(rows, columns=columns_ + ['extract_diagnosis', 'split_sent', 'histologic diagnosis'])

    # save null.csv
    result_df.loc[result_df['split_sent'].map(lambda x: True if x == '' else False)].to_csv(os.path.join(RESULT_PATH, 'null.csv'))
    result_df = result_df.loc[result_df['split_sent'].map(lambda x: True if x != '' else False)]
    result_df.reset_index(inplace=True, drop=True)
    result_df['first_'] = result_df['split_sent'].progress_apply(lambda x: ', '.join(x))

    tmp = result_df.copy()
    tmp['t'] = tmp['split_sent'].map(lambda x: [i for i in re.sub('(tissue labeled)|(specimen from)|(specimen labeled)|[\'\"]', '', ', '.join(x)).split(',')])
    tmp['labeling'] = tmp['t'].map(lambda x: labeling(x))
    tmp['organ'] = tmp['labeling'].map(lambda x: covert_list_to_str(x, 'organ'))
    tmp['location'] = tmp['labeling'].map(lambda x: covert_list_to_str(x, 'location'))
    tmp['opname'] = tmp['labeling'].map(lambda x: covert_list_to_str(x, 'opname'))
    tmp['histologic diagnosis'] = tmp['histologic diagnosis'].map(lambda x: re.sub('[->]', '', x).strip())
    tmp['direction'] = tmp['location'].map(lambda x: find_direction(x))
    tmp['lymph_node'] = tmp['histologic diagnosis'].map(lambda x: find_lymph_node(x))
    tmp['pre_lymph_node'] = tmp['lymph_node'].map(lambda x: preprocess_ln(x))

    df = tmp.loc[tmp['pre_lymph_node'].map(lambda x: True if len(x) % 2 == 0 else False)].copy()
    df.reset_index(inplace=True)

    lymphs = []

    for idx in range(df.shape[0]):
        lymph_info = []
        list_ = list()

        for i in df['pre_lymph_node'].iloc[idx]:

            if re.search('#', i):
                list_ = list()
                list_.append(i.strip())

            else:
                list_.append(i.strip())
                lymph_info.append(list_)

        lymphs.append(lymph_info)

    df['lymph_info'] = lymphs

    df_others = tmp.loc[tmp['pre_lymph_node'].map(lambda x: False if len(x) % 2 == 0 else True)].copy()
    df_others.reset_index(inplace=True)

    lymphs = []

    for idx in range(df_others.shape[0]):
        lymph_info = []
        list_ = list()

        for i in df_others['pre_lymph_node'].iloc[idx]:
            for j in i.split('#'):
                if len(j) >= 1:
                    if re.search('/', j):
                        list_.append(j.strip())
                        lymph_info.append(list_)

                    else:
                        if len(list_) == 1:
                            pre_ = list_[0]
                            list_ = list()
                            list_.append(pre_+' '+j)

                        else:
                            list_ = list()
                            if re.search('[0-9]', j):
                                list_.append('#'+j.strip())

        lymphs.append(lymph_info)

    df_others['lymph_info'] = lymphs

    df_total = pd.concat([df, df_others]).sort_values(by=['index']).set_index('index')
    df_total.reset_index(inplace=True, drop=True)

    lymphs = []

    for idx in range(df_total.shape[0]):
        lymph_info = []
        for i in df_total['lymph_info'].iloc[idx]:
            try:
                node_, positive_ = i[0].strip(), i[1]
                node_ = re.sub('(upper)|(lower)|(lobe)', '', node_)
                node_ = node_.replace('left', 'l')
                node_ = node_.replace('right', 'r')
                node_ = node_.replace('l&r', '')
                node_ = node_.replace('r&l', '')
                node_ = re.sub('[^#0-9rl&,\s]', '', node_)

                if (re.search('&', node_)) or (len(re.findall('#', node_)) > 1):
                    node_ = '#' + str(max([int(re.sub('[^0-9]', '', j)) for j in re.split('[&\s]', node_) if len(j) >= 1]))

                node_ = '#' + node_.split('#')[1]
                lymph_info.append([node_, positive_])

            except:
                pass

        lymphs.append(lymph_info)

    df_total['pre_lymph_info'] = lymphs
    df_total['tumor_size'] = df_total['histologic diagnosis'].map(lambda x: find_tumor_size(x))
    df_total['histologic diagnosis'] = df_total['histologic diagnosis'].map(lambda x: re.sub(',', ' ', re.split('[0-9\(\)]|(single)', x)[0]).strip())

    t_col = ['#{}'.format(i) for i in range(1, 14)] + ['#{}r'.format(i) for i in range(1, 14)] + ['#{}l'.format(i) for i in range(1, 14)]
    t = pd.DataFrame(columns=t_col, index=[i for i in range(df_total.shape[0])])

    for idx in range(df_total.shape[0]):
        for j in df_total['pre_lymph_info'].iloc[idx]:
            node_, positive_ = re.sub('\s', '', j[0]), j[1]

            try:
                if re.search('[rl]', node_):
                    t[node_].loc[idx]= positive_.replace('/', '|')
                else:
                    direction = list(set(df_total['direction'].iloc[idx]))
                    if len(direction) == 1:
                        node_ = node_ + direction[0][0]

                    t[node_].loc[idx]= positive_.replace('/', '|')

            except Exception:
                pass

    t = t.fillna('')

    total = pd.concat([df_total, t], axis=1)
    print('total shape: ', total.shape)

    # save total data
    total.to_csv(os.path.join(RESULT_PATH, 'final_result.csv'), index=False, encoding='utf-8-sig')

    # save final results
    total[columns_ +
           ['organ', 'location', 'opname',
           'histologic diagnosis', 'lymph_node', 'tumor_size',
           '#1', '#2', '#3', '#4', '#5', '#6', '#7', '#8', '#9', '#10', '#11', '#12', '#13',
           '#1r', '#2r','#3r', '#4r', '#5r', '#6r', '#7r', '#8r', '#9r', '#10r', '#11r', '#12r',
           '#13r', '#1l', '#2l', '#3l', '#4l', '#5l', '#6l', '#7l', '#8l', '#9l',
           '#10l', '#11l', '#12l', '#13l']].to_csv(os.path.join(RESULT_PATH, result_fname), index=False, encoding='utf-8-sig')