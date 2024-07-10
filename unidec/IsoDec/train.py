from IsoDec.engine import IsoDecEngine

if __name__ == "__main__":
    os.chdir("C:\\Data\\IsoNN\\multi")
    eng = IsoDecEngine(1)
    eng.create_training_dataloader("training_data_large.pth", "test_data_large.pth")
    eng.train_model(epochs=5)
