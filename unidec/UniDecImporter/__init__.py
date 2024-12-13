from unidec.UniDecImporter.ImporterFactory import ImporterFactory

def get_importer(file_path, **kwargs):
    return ImporterFactory.create_importer(file_path, **kwargs)

