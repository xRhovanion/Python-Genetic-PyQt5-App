# Python-Genetic-PyQt5-App

#### TÜRKÇE ####
PyQt5 tabanlı bir masaüstü uygulaması iskeleti: FASTA dosyası yükleme, sekans görüntüleme,
temel istatistikler (uzunluk, GC), nükleotid frekansı grafiği, motif arama (regex / sıradan arama / reverse-complement opsiyonları),
ORF bulucu + protein çevirisi ve restriksiyon enzimleri (EcoRI, HindIII vb.) kesim bölgelerini tespit etme özellikleri içerir.

Amaç:
Kapsamlı bir genetik analiz uygulamasının ilk modüllerini hızlıca prototiplemek; bu dosya
üyelikli modüller (motif arama, ORF bulucu, CRISPR analizi, görselleştirme, rapor) ile genişletilebilir.

Gereksinimler:

Python 3.8+

PyQt5

matplotlib

numpy

(opsiyonel) biopython — varsa FASTA parsing için kullanılır, yoksa internal parser çalışır.

Kurulum örneği:

pip install PyQt5 matplotlib numpy biopython qdarkstyle


Çalıştırma:

python genetik_pyqt5_app.py

#### ENGLISH ####
A PyQt5-based desktop application skeleton: load FASTA files, view sequences,
basic statistics (length, GC content), nucleotide frequency chart, motif search (regex / simple search / reverse-complement options),
ORF finder + protein translation, and restriction enzyme (e.g., EcoRI, HindIII) cut site detection.

Purpose:
To quickly prototype the initial modules of a comprehensive genetic analysis application;
this file can be extended with modular add-ons (motif search, ORF finder, CRISPR analysis, visualization, reporting).

Requirements:

Python 3.8+

PyQt5

matplotlib

numpy

(optional) biopython — used for FASTA parsing if available, otherwise internal parser runs.

Installation example:
    
pip install PyQt5 matplotlib numpy biopython qdarkstyle


Run:

python genetik_pyqt5_app.py
