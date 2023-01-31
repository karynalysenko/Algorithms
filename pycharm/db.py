import mysql.connector as SQLC
class DB:
    # database connection and extraction of sequences and id_genebanks
    def __init__(self):
        self.cnx = SQLC.connect(
            host="geo.di.uminho.pt",
            user="bioinformatica",
            password="20221207",
            database="AP_db_KRG")
        self.sequences = []

    def extract_sequences(self):
        cursor = self.cnx.cursor()
        cursor.execute("SELECT sequence FROM Gene")
        self.sequences = cursor.fetchall()
        cursor.close()

    def extract_ids(self):
        cursor = self.cnx.cursor()
        cursor.execute("SELECT ID_genebank FROM Gene")
        self.sequences = cursor.fetchall()
        cursor.close()

    def __str__(self):
        return self.sequences