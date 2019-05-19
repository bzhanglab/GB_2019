
from numpy.random import seed
seed(2019)
from tensorflow import set_random_seed
set_random_seed(2020)

import keras
from keras.callbacks import EarlyStopping, ModelCheckpoint, CSVLogger, ReduceLROnPlateau
from keras.layers import Dense, Dropout, Activation, Flatten, Input, MaxPooling2D, Conv2D, Conv1D, Bidirectional, LSTM, \
    Embedding, MaxPooling1D, Average, CuDNNGRU, CuDNNLSTM, Bidirectional, GRU
from keras.layers.normalization import BatchNormalization
from keras.models import Model, load_model, clone_model, model_from_json
from keras.models import Sequential
from keras.optimizers import Adam, SGD

from sklearn.model_selection import train_test_split

import os
import scipy
from keras.callbacks import Callback
import sklearn

import pickle

import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns

import argparse
import sys

from math import ceil, floor

import os
import gc
import keras.backend as K
import tensorflow as tf
import time
import datetime

import json
from shutil import copyfile





letterDict = {"A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7, "K": 8, "L": 9, "M": 10, "N": 11,
                  "P": 12, "Q": 13, "R": 14, "S": 15, "T": 16, "V": 17, "W": 18, "Y": 19, "U": 20, "B": 21, "X": 22}

def add_mod(mod=None):

    i=len(letterDict)
    for m in mod:
        if str(m) not in letterDict:
            letterDict[str(m)] = i
            print("New aa: %s -> %d" % (str(m),i))
            i = i + 1

def load_aa(file):
    print("Load aa coding data from file %s" % (file))
    dat = pd.read_table(file, sep="\t", header=0, low_memory=False)
    letterDict.clear()
    for i, row in dat.iterrows():
        letterDict[row['aa']] = row['i']

def save_aa(file):
    print("Save aa coding data to file %s" % (file))
    with open(file,"w") as f:
        f.write("aa\ti\n")
        for aa in letterDict.keys():
            f.write(aa+"\t"+str(letterDict[aa])+"\n")


def data_processing(input_data: str, test_file=None, mod=None, max_x_length = 50, min_rt=0, max_rt=120, unit="s",
                    out_dir="./", aa_file=None):

    if aa_file is not None:
        ## read aa information from file
        load_aa(aa_file)
    else:
        if mod is not None:
            add_mod(mod)

        aa2file = out_dir + "/aa.tsv"
        save_aa(aa2file)

    ##
    siteData = pd.read_table(input_data, sep="\t", header=0, low_memory=False)

    if "x" not in siteData.columns:
        siteData.columns = ['x','y']

    if unit.startswith("s"):
        siteData['y'] = siteData['y']/60.0

    if max_rt < siteData['y'].max():
        max_rt = siteData['y'].max() + 1.0

    if min_rt > siteData['y'].min():
        min_rt = siteData['y'].min() - 1.0

    # aaMap = getAAcodingMap()
    n_aa_types = len(letterDict)
    print("AA types: %d" % (n_aa_types))


    ## all aa in data
    all_aa = set()

    ## get the max length of input sequences

    longest_pep_training_data = 0
    for pep in siteData["x"]:
        if max_x_length < len(pep):
            max_x_length = len(pep)

        if longest_pep_training_data < len(pep):
            longest_pep_training_data = len(pep)

        ##
        for aa in pep:
            all_aa.add(aa)



    print("Longest peptide in training data: %d\n" % (longest_pep_training_data))

    ## test data
    test_data = None
    longest_pep_test_data = 0
    if test_file is not None:
        print("Use test file %s" % (test_file))
        test_data = pd.read_table(test_file, sep="\t", header=0, low_memory=False)
        if "x" not in test_data.columns:
            test_data.columns = ['x', 'y']
        if unit.startswith("s"):
            test_data['y'] = test_data['y'] / 60.0

        if max_rt < test_data['y'].max():
            max_rt = test_data['y'].max() + 1.0

        if min_rt > test_data['y'].min():
            min_rt = test_data['y'].min() - 1.0

        for pep in test_data["x"]:
            if max_x_length < len(pep):
                max_x_length = len(pep)

            if longest_pep_test_data < len(pep):
                longest_pep_test_data = len(pep)

            for aa in pep:
                all_aa.add(aa)

        print("Longest peptide in test data: %d\n" % (longest_pep_test_data))

    print(sorted(all_aa))

    siteData = siteData.sample(siteData.shape[0], replace=False, random_state=2018)

    train_data = np.zeros((siteData.shape[0], max_x_length, n_aa_types))
    k = 0
    for i, row in siteData.iterrows():
        peptide = row['x']
        train_data[k] = encodePeptideOneHot(peptide, max_length=max_x_length)
        k = k + 1

    train_data = train_data.reshape(train_data.shape[0], train_data.shape[1], train_data.shape[2])

    X_test = np.empty(1)
    Y_test = np.empty(1)

    print("RT range: %d - %d\n" % (min_rt,max_rt))

    if test_data is None:
        X_train, X_test, Y_train, Y_test = train_test_split(train_data,
                                                            #to_categorical(pos_neg_all_data['y'], num_classes=2),
                                                            minMaxScale(siteData['y'],min_rt,max_rt),
                                                            test_size=0.1, random_state=100)
    else:
        X_train = train_data
        #Y_train = to_categorical(pos_neg_all_data['y'], num_classes=2)
        Y_train = siteData['y']
        Y_train = minMaxScale(Y_train, min_rt, max_rt)
        if len(Y_train.shape) >= 2:
            Y_train = Y_train.reshape(Y_train.shape[1])

        X_test = np.zeros((test_data.shape[0], max_x_length, n_aa_types))
        k = 0
        for i, row in test_data.iterrows():
            peptide = row['x']
            X_test[k] = encodePeptideOneHot(peptide, max_length=max_x_length)
            k = k + 1

        Y_test = minMaxScale(test_data['y'],min_rt,max_rt)


    X_train = X_train.astype('float32')
    X_test = X_test.astype('float32')

    print("X_train shape:")
    print(X_train.shape)
    print("X_test shape:")
    print(X_test.shape)
    print("Modeling start ...")


    return [X_train, Y_train, X_test, Y_test, min_rt, max_rt]


def processing_prediction_data(model_file: str, input_data: str):
    '''

    :param model_file: model file in json format
    :param input_data: prediction file
    :return: A numpy matrix for prediction
    '''

    with open(model_file, "r") as read_file:
        model_list = json.load(read_file)

    model_folder = os.path.dirname(model_file)
    aa_file = model_folder + "/" + os.path.basename(model_list['aa'])
    load_aa(aa_file)

    ##
    siteData = pd.read_table(input_data, sep="\t", header=0, low_memory=False)

    n_aa_types = len(letterDict)
    print("AA types: %d" % (n_aa_types))

    ## all aa in data
    all_aa = set()

    ## get the max length of input sequences
    max_x_length = model_list['max_x_length']

    longest_pep_len = 0
    for pep in siteData["x"]:
        #if max_x_length < len(pep):
        #    max_x_length = len(pep)
        if longest_pep_len < len(pep):
            longest_pep_len = len(pep)
        ##
        for aa in pep:
            all_aa.add(aa)

    print("Longest peptide in input data: %d\n" % (longest_pep_len))

    print(sorted(all_aa))

    # siteData = siteData.sample(siteData.shape[0], replace=False, random_state=2018)

    train_data = np.zeros((siteData.shape[0], max_x_length, n_aa_types))
    k = 0
    for i, row in siteData.iterrows():
        peptide = row['x']
        train_data[k] = encodePeptideOneHot(peptide, max_length=max_x_length)
        k = k + 1

    train_data = train_data.reshape(train_data.shape[0], train_data.shape[1], train_data.shape[2])


    train_data = train_data.astype('float32')

    return train_data



def build_default_model(input_shape):
    model = Sequential()
    model.add(Conv1D(512, 3, padding='same', input_shape=input_shape))
    model.add(Activation('relu'))
    model.add(BatchNormalization())
    model.add(MaxPooling1D(pool_size=(2)))

    model.add(Conv1D(128, 3, padding='same'))
    model.add(Activation('relu'))
    model.add(BatchNormalization())
    model.add(MaxPooling1D(pool_size=(2)))

    model.add(Conv1D(64, 3, padding='same'))
    model.add(Activation('relu'))
    model.add(BatchNormalization())
    model.add(MaxPooling1D(pool_size=(2)))

    # model.add(Bidirectional(CuDNNGRU(50, return_sequences=True)))
    model.add(Bidirectional(GRU(50, return_sequences=True)))
    # model.add(SeqSelfAttention(attention_activation='sigmoid'))
    model.add(Dropout(0.5))

    model.add(Flatten())

    model.add(Dense(64))
    model.add(Activation('relu'))
    model.add(Dropout(0.45))

    model.add(Dense(1, kernel_initializer='normal', activation='sigmoid'))
    return model

def train_tf(input_data: str, test_file=None, batch_size=64, nb_epoch=100, early_stop=None, mod=None,
                max_x_length = 50, min_rt=0, max_rt=120, unit="s",out_dir="./", prefix = "test",
                model_file=None):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    print("Build deep learning model ...")

    X_train, Y_train, X_test, Y_test, min_rt, max_rt = data_processing(input_data=input_data, test_file = test_file, mod = mod, max_x_length = max_x_length,
                    min_rt = min_rt, max_rt = max_rt, unit = unit)





def train_model(input_data: str, test_file=None, batch_size=64, nb_epoch=100, early_stop=None, mod=None,
                max_x_length = 50, min_rt=0, max_rt=120, unit="s",out_dir="./", prefix = "test",
                p_model=None,
                model=None, optimizer_name=None):

    res_map = dict()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    print("Build deep learning model ...")

    X_train, Y_train, X_test, Y_test, min_rt, max_rt = data_processing(input_data=input_data, test_file = test_file, mod = mod, max_x_length = max_x_length,
                    min_rt = min_rt, max_rt = max_rt, unit = unit)

    if model is None:
        print("Use default model ...")
        model = build_default_model(X_train.shape[1:])
    else:
        print("Use input model ...")
        model = clone_model(model)

    if p_model is not None:
        transfer_layer = 5
        frozen = True
        # model_copy.set_weights(model.get_weights())
        base_model = load_model(p_model)
        print("Perform transfer learning ...")
        n_layers = len(base_model.layers)
        print("The number of layers: %d" % (n_layers))
        for l in range((n_layers - transfer_layer)):
            if l != 0:
                model.layers[l].set_weights(base_model.layers[l].get_weights())
                if frozen is True:
                    model.layers[l].trainable = False
                    print("layer (frozen:True): %d %s" % (l,model.layers[l].name))
                else:
                    print("layer (frozen:False): %d %s" % (l,model.layers[l].name))

    if model.optimizer is None:
        ## use default optimizer: Adam
        if optimizer_name is None:
            print("Use default optimizer:Adam")
            model.compile(loss='mean_squared_error',
                          optimizer="adam",
                          metrics=['mse', 'mae'])
        else:
            print("Use optimizer provided by user: %s" % (optimizer_name))
            model.compile(loss='mean_squared_error',
                          optimizer=optimizer_name,
                          metrics=['mse', 'mae'])

    else:
        if optimizer_name is None:
            print("Use optimizer from the model.")
            model.compile(loss='mean_squared_error',
                               ## In this case, we cannot change the learning rate.
                               optimizer=model.optimizer,
                               metrics=['mse','mae'])


        else:
            print("Use optimizer provided by user: %s" % (optimizer_name))
            model.compile(loss='mean_squared_error',
                               ## In this case, we cannot change the learning rate.
                               optimizer=optimizer_name,
                               metrics=['mse','mae'])


    print("optimizer: %s" % (type(model.optimizer)))

    model.summary()
    # model = multi_gpu_model(model, gpus=3)

    my_callbacks = RegCallback(X_train, X_test, Y_train, Y_test, min_rt=min_rt, max_rt=max_rt)
    # Save model
    model_chk_path = out_dir + "/best_model.hdf5"
    mcp = ModelCheckpoint(model_chk_path, monitor="val_mean_squared_error", save_best_only=True, save_weights_only=False,
                          verbose=1, mode='min')

    ## monitor training information
    # tbCallBack = callbacks.TensorBoard(log_dir='./Graph', histogram_freq=0, write_graph=True, write_images=True)
    #model.fit(X_train, Y_train, batch_size=batch_size, epochs=nb_epoch, validation_data=(X_test, Y_test), callbacks=[my_callbacks, mcp])
    model.fit(X_train, Y_train, batch_size=batch_size, epochs=nb_epoch, validation_split=0.1,
              callbacks=[my_callbacks, mcp])

    ## get the best model
    model_best = load_model(model_chk_path)

    y_pred = model_best.predict(X_test)

    y_pred_rev = minMaxScoreRev(y_pred, min_rt, max_rt)
    y_true = minMaxScoreRev(Y_test, min_rt, max_rt)

    #test_data['pred'] = y_pred_rev
    #test_data.to_csv("pred.csv")
    x = pd.DataFrame({"y": y_true, "y_pred": y_pred_rev.reshape(y_pred_rev.shape[0])})
    out_file = out_dir + "/" + prefix +".csv"
    print("Prediction result: %s" % (out_file))
    x.to_csv(out_file)

    res_map['model'] = model_best
    return res_map



def minMaxScale(x, min=0,max=120):
    new_x = 1.0*(x-min)/(max-min)
    return new_x

def minMaxScoreRev(x,min=0,max=120):
    old_x = x * (max - min) + min
    return old_x



def encodePeptideOneHot(peptide: str, max_length=None):  # changed add one column for '1'

    AACategoryLen = len(letterDict)
    peptide_length = len(peptide)
    use_peptide = peptide
    if max_length is not None:
        if peptide_length < max_length:
            use_peptide = peptide + "X" * (max_length - peptide_length)

    en_vector = np.zeros((len(use_peptide), AACategoryLen))

    i = 0
    for AA in use_peptide:
        if AA in letterDict.keys():
            try:
                en_vector[i][letterDict[AA]] = 1
            except:
                print("peptide: %s, i => aa: %d, %s, %d" % (use_peptide,i, AA, letterDict[AA]))
                exit(1)
        else:
            en_vector[i] = np.full(AACategoryLen,1/AACategoryLen)

        i = i + 1

    return en_vector

def ensemble_models(input_data: str, test_file=None,models=None,
                    models_file=None,
                    genome_handler_file=None,
                    top_n=10,
                    trained=True,
                    ensemble_method="average",
                    batch_size=64, nb_epoch=100, early_stop=None, mod=None,
                    max_x_length=50, min_rt=0, max_rt=120, unit="s", out_dir="./", prefix="test"):
    """
    This function is used to ensemble multiple deep learning models. It can be used for training and testing.
    For each model, predict
    1. testing
    2. prediction
    3. training and testing:
    a=ensemble_models(models_file="genomes.csv",genome_handler_file="genome_handler.obj",top_n=4,input_file="filtered_input_data.txt",test_file="final_train.txt")
    4. training
    :return:
    """
    from AutoSeq import GenomeHandler

    # print("The number of models:", len(models))

    # test data
    X_test = np.empty(1)
    Y_test = np.empty(1)

    y_pr = []
    score = []

    model_list = dict()


    if genome_handler_file is not None:
        X_train, Y_train, X_test, Y_test, min_rt, max_rt = data_processing(input_data=input_data, test_file=test_file,
                                                                           mod=mod, max_x_length=max_x_length,
                                                                           min_rt=min_rt, max_rt=max_rt, unit=unit,
                                                                           out_dir=out_dir)
        model_list['dp_model'] = dict()
        model_list['max_x_length'] = X_train.shape[1]
        model_list['aa'] = out_dir + "/aa.tsv"
        print("max_x_length: %s" % (max_x_length))
        # read models from genetic search result configure file
        optimizer_name = dict()
        if models_file is not None:
            models = dict()
            gn = pd.read_csv(models_file)
            select_models = gn.sort_values('Val Accuracy', ascending=True).head(top_n)
            genome_handler = pickle.load(open(genome_handler_file, "rb"))
            genome_handler.input_shape = X_train.shape[1:]
            select_models = np.array(select_models.iloc[:, 0:(select_models.shape[1] - 2)])
            for i in range(0, select_models.shape[0]):
                #models[i], optimizer_name = genome_handler.decodeOneHot(select_models[i],return_optimizer=True)
                models[i], optimizer_name[i] = genome_handler.decodeOneHotPlusLSTM(select_models[i], return_optimizer=True)

            trained = False
        else:
            print("")

        if not trained:
            print("Training ...")
            # For each model, train the model
            for (name, model) in models.items():
                print("Train model:", name)
                # perform sample specific training
                res_map = train_model(input_data=input_data, test_file=test_file, batch_size=batch_size,
                                      nb_epoch=nb_epoch, early_stop=early_stop, mod=mod,
                                      max_x_length=max_x_length, min_rt=min_rt, max_rt=max_rt, unit=unit,
                                      out_dir=out_dir, prefix=str(name), model=model,
                                      optimizer_name=optimizer_name[name])

                ## save the model to a file:
                model_file_name = "model_" + str(name) + ".h5"
                model_file_path = out_dir + "/" + model_file_name
                res_map["model"].save(model_file_path)

                model_list['dp_model'][name] = model_file_path

                del res_map
                gc.collect()
                K.clear_session()
                tf.reset_default_graph()
        else:
            print("The models have been trained!")


    else:

        ## Transfer learning
        with open(models_file, "r") as read_file:
            model_list = json.load(read_file)

        model_folder = os.path.dirname(models_file)
        aa_file = os.path.basename(model_list['aa'])
        aa_file = model_folder + "/" + aa_file
        X_train, Y_train, X_test, Y_test, min_rt, max_rt = data_processing(input_data=input_data, test_file=test_file,
                                                                           mod=mod, max_x_length=model_list['max_x_length'],
                                                                           min_rt=min_rt, max_rt=max_rt, unit=unit,
                                                                           out_dir=out_dir,aa_file=aa_file)


        new_model_list = dict()
        new_model_list['dp_model'] = dict()
        for (name, dp_model_file) in model_list['dp_model'].items():
            print("\nDeep learning model:", name)
            # keras model evaluation: loss and accuracy
            # load model
            model_name = os.path.basename(dp_model_file)
            model_full_path = model_folder + "/" + model_name

            model = load_model(model_full_path)
            #new_model = change_model(model, X_train.shape[1:])
            new_model = model

            print("Perform transfer learning ...")
            n_layers = len(new_model.layers)
            print("The number of layers: %d" % (n_layers))
            #for layer in new_model.layers:
            #    layer_name = str(layer.name)
            #    if layer_name.startswith("dense"):
            #        break
            #    else:
            #        layer.trainable = False
            #        print("layer (frozen:True): %s" % (layer_name))

            new_model.compile(loss='mean_squared_error',
                              ## In this case, we cannot change the learning rate.
                              optimizer=model.optimizer,
                              #optimizer=Adam(lr=0.0001),
                              #optimizer=SGD(lr=1e-3, decay=1e-4, momentum=0.9, nesterov=True),
                              metrics=['mse', 'mae'])
            my_callbacks = RegCallback(X_train, X_test, Y_train, Y_test, min_rt=min_rt, max_rt=max_rt)
            # Save model
            model_chk_path = out_dir + "/best_model.hdf5"
            mcp = ModelCheckpoint(model_chk_path, monitor="val_mean_squared_error", save_best_only=True,
                                  save_weights_only=False,
                                  verbose=1, mode='min')

            ## monitor training information
            # tbCallBack = callbacks.TensorBoard(log_dir='./Graph', histogram_freq=0, write_graph=True, write_images=True)
            new_model.fit(X_train, Y_train, batch_size=batch_size, epochs=nb_epoch, validation_data=(X_test, Y_test),
                      callbacks=[my_callbacks, mcp])

            ## get the best model
            best_model = load_model(model_chk_path)
            ## save the model to a file:
            model_file_name = "model_" + str(name) + ".h5"
            model_file_path = out_dir + "/" + model_file_name
            best_model.save(model_file_path)

            new_model_list['dp_model'][name] = model_file_path

            gc.collect()
            K.clear_session()
            tf.reset_default_graph()

        new_model_list['max_x_length'] = model_list['max_x_length']
        new_aa_file = out_dir + "/" + os.path.basename(model_list['aa'])
        copyfile(aa_file, new_aa_file)
        new_model_list['aa'] = new_aa_file

        ## Useful for new data prediction
        new_model_list['min_rt'] = min_rt
        new_model_list['max_rt'] = max_rt

        model_list = new_model_list


    # save model data
    #file_all_models = open(out_dir + "/all_models.obj", 'wb')
    #pickle.dump(models, file_all_models)
    #file_all_models.close()

    ####################################################################################################################
    print("Ensemble learning ...")


    para = dict()
    para['min_rt'] = min_rt
    para['max_rt'] = max_rt

    ## save result
    model_json = out_dir + "/model.json"
    with open(model_json, 'w') as f:
        json.dump(model_list, f)

    ## evaluation
    if test_file is not None:
        ensemble_predict(model_json,x=X_test,y=Y_test,para=para, batch_size=batch_size,method=ensemble_method,
                         out_dir=out_dir,
                         prefix="final_eval")

    ####################################################################################################################

def change_model(model, new_input_shape):
    # replace input shape of first layer
    print("Base model ...")
    print(model.get_weights())
    model._layers[1].batch_input_shape = (None,new_input_shape[0],new_input_shape[1])

    # rebuild model architecture by exporting and importing via json
    new_model = keras.models.model_from_json(model.to_json())


    # copy weights from old model to new one
    for layer in new_model.layers:
        try:
            print("layer: %s" % (layer.name))
            layer.set_weights(model.get_layer(name=layer.name).get_weights())
        except:
            print("Could not transfer weights for layer {}".format(layer.name))

    new_model.compile(loss='mean_squared_error',
                  ## In this case, we cannot change the learning rate.
                  optimizer=model.optimizer,
                  metrics=['mse', 'mae'])

    new_model.summary()

    print("New model ...")
    print(new_model.get_weights())
    return new_model


def ensemble_predict(model_file:str, x, y=None, para=None, out_dir="./", method="average",batch_size=64, prefix="test"):

    res_to_file = np.empty(0)
    ## prediction result
    y_pr_final = np.empty(0)

    if method == "average":
        print("Average ...")
        res = dl_models_predict(model_file, x=x, y=y, para=para, batch_size=batch_size,out_dir=out_dir,prefix=prefix)
        y_pr_final = res.mean(axis=1)
        #np.save("res", res)
        #np.save("y_pr_final",y_pr_final)
        res_to_file = res
        res_to_file = np.append(res_to_file, y_pr_final.reshape([y_pr_final.shape[0],1]), axis=1)

    if y is not None:
        # evaluate the final _result
        # ROC
        print("\n\nFinal model:")
        if len(y.shape) >= 2:
            y_true_class = np.argmax(y, axis=1)
        else:
            y_true_class = y

        # classification report: precision recall f1-score support
        # y_class_final = np.where(y_pr_final > 0.5, 1, 0) #np.argmax(y_pr_final, axis=1)
        out_prefix = prefix + "_final"
        evaluate_model(y_true_class, y_pr_final, para=para, plot=True, out_dir=out_dir, prefix=out_prefix)

        #np.save("y",y)
        #res_to_file = np.append(res_to_file, y_true_class.reshape([y_true_class.shape[0],1]), axis=1)

    #np.save("final_res",res_to_file)
    if y is not None:
        return res_to_file
    else:
        return y_pr_final

def rt_predict(model_file:str, test_file:str, out_dir="./", prefix="test", method = "average"):

    res_to_file = np.empty(0)
    ## prediction result
    y_pr_final = np.empty(0)

    x_predict_data = processing_prediction_data(model_file, test_file)

    if method == "average":
        print("Average ...")
        res = dl_models_predict(model_file, x=x_predict_data, out_dir=out_dir, prefix=prefix)
        y_pr_final = res.mean(axis=1)
        #np.save("res", res)
        #np.save("y_pr_final", y_pr_final)
        res_to_file = res
        #res_to_file = np.append(res_to_file, y_pr_final.reshape([y_pr_final.shape[0], 1]), axis=1)

        with open(model_file, "r") as read_file:
            model_list = json.load(read_file)

        rt_pred = minMaxScoreRev(y_pr_final, model_list['min_rt'], model_list['max_rt'])

        input_data = pd.read_table(test_file, sep="\t", header=0, low_memory=False)
        input_data['y_pred'] = rt_pred

        ## output
        out_file = out_dir + "/" + prefix + ".csv"
        input_data.to_csv(out_file,sep="\t",index=False)





def dl_models_predict(model_file, x, y=None, para=None,batch_size=64, out_dir="./", prefix="test"):

    with open(model_file, "r") as read_file:
        model_list = json.load(read_file)


    y_dp = np.zeros(0)

    model_folder = os.path.dirname(model_file)
    avg_models = list()
    for (name, dp_model_file) in model_list['dp_model'].items():
        print("\nDeep learning model:", name)
        # keras model evaluation: loss and accuracy
        # load model
        model_name = os.path.basename(dp_model_file)
        model_full_path = model_folder + "/" + model_name

        model = load_model(model_full_path)
        avg_models.append(model)
        y_prob = model.predict(x, batch_size=batch_size)
        ## for class 1
        #y_prob_dp_vector = y_prob[:, 1]
        y_prob_dp_vector = y_prob
        y_prob_dp = y_prob_dp_vector.reshape([y_prob_dp_vector.shape[0], 1])
        if y_dp.shape[0] != 0:
            y_dp = np.append(y_dp, y_prob_dp, axis=1)
        else:
            y_dp = y_prob_dp

        if y is not None:
            evaluation_res = model.evaluate(x, y)
            print("Metrics:")
            print(evaluation_res)

            # ROC
            if len(y.shape) >= 2:
                y_true_class = np.argmax(y, axis=1)
            else:
                y_true_class= y
            out_prefix = prefix + "_" + str(name)
            evaluate_model(y_true_class, y_prob_dp_vector, para=para, plot=True, out_dir=out_dir, prefix=out_prefix)


        gc.collect()
        K.clear_session()
        tf.reset_default_graph()

    return y_dp




def evaluate_model(y_t, y_p, para=None, plot=True, out_dir="./", prefix="test"):
    y_t = minMaxScoreRev(y_t, para['min_rt'], para['max_rt'])
    y_p = minMaxScoreRev(y_p, para['min_rt'], para['max_rt'])
    y2 = pd.DataFrame({"y": y_t, "y_pred": y_p.reshape(y_p.shape[0])})
    cor = scipy.stats.pearsonr(y2['y'], y2['y_pred'])[0]
    mae = sklearn.metrics.mean_absolute_error(y2['y'], y2['y_pred'])
    r2 = sklearn.metrics.r2_score(y2['y'], y2['y_pred'])
    abs_median = float(np.median(np.abs(y2['y'] - y2['y_pred'])))
    d_t95 = Delta_t95(y2['y'], y2['y_pred'])
    print('Cor: %s, MAE: %s, R2: %s, abs_median_e: %s, dt95: %s' % (
        str(round(cor, 4)), str(round(mae, 4)),
        str(round(r2, 4)), str(round(abs_median, 4)),
        str(round(d_t95, 4))), end=100 * ' ' + '\n')

    ## output
    out_file = out_dir + "/" + prefix + ".csv"
    y2.to_csv(out_file)


class RegCallback(Callback):
    """
    Calculate AUROC for each epoch
    """

    def __init__(self, X_train, X_test, y_train, y_test, min_rt=0, max_rt=120):
        self.x = X_train
        self.y = minMaxScoreRev(y_train,min_rt,max_rt)
        self.x_val = X_test
        self.y_val = minMaxScoreRev(y_test,min_rt,max_rt)
        self.min_rt = min_rt
        self.max_rt = max_rt

    def on_train_begin(self, logs={}):
        return

    def on_train_end(self, logs={}):
        return

    def on_epoch_begin(self, epoch, logs={}):
        return

    def on_epoch_end(self, epoch, logs=None):

        ## training data
        y_pred = self.model.predict(self.x)
        y_pred_rev = minMaxScoreRev(y_pred, self.min_rt, self.max_rt)

        y1 = pd.DataFrame({"y": self.y, "y_pred": y_pred_rev.reshape(y_pred_rev.shape[0])})
        cor1 = scipy.stats.pearsonr(y1['y'],y1['y_pred'])[0]
        mae1 = sklearn.metrics.mean_absolute_error(y1['y'],y1['y_pred'])
        r21 = sklearn.metrics.r2_score(y1['y'],y1['y_pred'])
        abs_median1 = np.median(np.abs(y1['y'] - y1['y_pred']))
        d_t951 = Delta_t95(y1['y'], y1['y_pred'])
        ## test data
        y_pred_val = self.model.predict(self.x_val)
        y_pred_val_rev = minMaxScoreRev(y_pred_val, self.min_rt, self.max_rt)
        y2 = pd.DataFrame({"y": self.y_val, "y_pred": y_pred_val_rev.reshape(y_pred_val_rev.shape[0])})
        cor2 = scipy.stats.pearsonr(y2['y'], y2['y_pred'])[0]
        mae2 = sklearn.metrics.mean_absolute_error(y2['y'], y2['y_pred'])
        r22 = sklearn.metrics.r2_score(y2['y'], y2['y_pred'])
        abs_median2 = np.median(np.abs(y2['y'] - y2['y_pred']))
        d_t952 = Delta_t95(y2['y'], y2['y_pred'])
        print('\nCor: %s - Cor_val: %s, MAE: %s - MAE_val: %s, R2: %s - R2_val: %s, MedianE: %s - MedianE_val: %s, dt95: %s - dt95_val: %s' % (str(round(cor1, 4)), str(round(cor2, 4)),
                                                                                       str(round(mae1, 4)), str(round(mae2, 4)),
                                                                                       str(round(r21, 4)), str(round(r22, 4)),
                                                                                       str(round(abs_median1, 4)), str(round(abs_median2, 4)),
                                                                                       str(round(d_t951, 4)), str(round(d_t952, 4))), end=100 * ' ' + '\n')

        return

    def on_batch_begin(self, batch, logs={}):
        return

    def on_batch_end(self, batch, logs={}):
        return

def Delta_t95(obs, pred):
    q95 = int(np.ceil(len(obs) * 0.95))
    return 2 * sorted(abs(obs - pred))[q95 - 1]


def main():

    if len(sys.argv) == 1:
        print("python autort.py [train, predict]")
        sys.exit(0)
    else:

        mode = sys.argv[1]

        if mode == "train":

            ## general training and training for transfer learning
            parser = argparse.ArgumentParser(
                description='AutoRT')
            parser.add_argument('-i', '--input', default=None, type=str, required=True,
                                help="Input data for training")
            parser.add_argument('-t', '--test', default=None, type=str,
                                help="Input data for testing")

            parser.add_argument('-o', '--out_dir', default="./", type=str,
                                help="Output directory")

            parser.add_argument('-e', '--epochs', default=20, type=int)
            parser.add_argument('-b', '--batch_size', default=128, type=int)
            parser.add_argument('-r2', '--max_rt', default=0, type=int)
            parser.add_argument('-l', '--max_length', default=0, type=int)
            parser.add_argument('-m', '--mod', default=None, type=str)
            parser.add_argument('-u', '--unit', default="s", type=str)

            parser.add_argument('-s', '--ensemble', default=None, type=str)
            parser.add_argument('-s2', '--ga', default=None, type=str)
            parser.add_argument('-w', '--top_n', default=10, type=int)

            args = parser.parse_args(sys.argv[2:len(sys.argv)])

            input_file = args.input
            test_file = args.test
            out_dir = args.out_dir
            max_rt = args.max_rt
            max_length = args.max_length
            mod = args.mod
            unit = args.unit

            if mod is not None:
                mod = mod.split(",")

            epochs = args.epochs
            batch_size = args.batch_size

            ensemble = args.ensemble
            top_n = args.top_n
            ga = args.ga

            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            ensemble_models(input_data=input_file, test_file=test_file, nb_epoch=epochs, batch_size=batch_size,
                                top_n=top_n,
                                max_rt=max_rt, max_x_length=max_length, mod=mod, unit=unit, models_file=ensemble,
                                genome_handler_file=ga,
                                out_dir=out_dir)

        elif mode == "predict":

            parser = argparse.ArgumentParser(
                description='AutoRT')

            parser.add_argument('-t', '--test', default=None, type=str,
                                help="Input data for testing")
            parser.add_argument('-o', '--out_dir', default="./", type=str,
                                help="Output directory")
            parser.add_argument('-p', '--prefix', default="test", type=str)
            parser.add_argument('-s', '--ensemble', default=None, type=str)

            args = parser.parse_args(sys.argv[2:len(sys.argv)])

            test_file = args.test
            out_dir = args.out_dir
            prefix = args.prefix


            ensemble = args.ensemble

            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            rt_predict(model_file=ensemble,test_file=test_file, out_dir=out_dir, prefix=prefix)




if __name__=="__main__":
    main()






