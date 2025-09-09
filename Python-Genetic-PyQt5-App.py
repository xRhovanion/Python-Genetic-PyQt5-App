"""
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

Python-Genetic-PyQt5-App.py

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

Python-Genetic-PyQt5-App.py
"""

import sys
import os
import re
import csv
from collections import Counter

try:
    from PyQt5.QtWidgets import (
        QApplication, QMainWindow, QFileDialog, QListWidget, QTextEdit, QTabWidget,
        QWidget, QVBoxLayout, QPushButton, QLabel, QLineEdit, QHBoxLayout, QAction,
        QSplitter, QTableWidget, QTableWidgetItem, QMessageBox, QCheckBox, QSpinBox,
        QComboBox
    )
    from PyQt5.QtCore import Qt
    from PyQt5.QtGui import QTextCursor, QTextCharFormat, QColor
except Exception as e:
    raise RuntimeError("PyQt5 gerekli. `pip install PyQt5` ile yükleyin. Hata: {}".format(e))

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Biopython
try:
    from Bio import SeqIO
    BIOPYTHON = True
except Exception:
    BIOPYTHON = False


# ---------------------- Helper functions ----------------------

def parse_fasta(path):
    """Basit FASTA parser. Eğer Biopython varsa onu kullanır."""
    """Simple FASTA parser. Uses Biopython if available."""
    sequences = {}
    if BIOPYTHON:
        for rec in SeqIO.parse(path, "fasta"):
            sequences[str(rec.id)] = str(rec.seq).upper().replace('\n', '').replace(' ', '')
        return sequences

    with open(path, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    sequences[header] = ''.join(seq_lines).upper()
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            sequences[header] = ''.join(seq_lines).upper()
    return sequences


def revcomp(seq):
    trans = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(trans)[::-1]


def gc_content(seq):
    if len(seq) == 0:
        return 0.0
    seq = seq.upper()
    g = seq.count('G')
    c = seq.count('C')
    return (g + c) / len(seq) * 100


def nucleotide_freq(seq):
    seq = seq.upper()
    counts = Counter(seq)
    total = sum(counts[nuc] for nuc in ['A', 'C', 'G', 'T'])
    freqs = {nuc: counts.get(nuc, 0) for nuc in ['A', 'C', 'G', 'T']}
    return freqs, total


def sliding_gc(seq, window=50):
    seq = seq.upper()
    L = len(seq)
    if L == 0 or window <= 0:
        return [], []
    gcs = []
    positions = []
    for i in range(0, L - window + 1):
        win = seq[i:i+window]
        gcs.append(gc_content(win))
        positions.append(i + window//2)
    return positions, gcs


# Minimal genetic code
CODON_TABLE = {
    # Standard codon table (one-letter aa)
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}


def translate_dna(seq):
    seq = seq.upper()
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        protein.append(aa)
    return ''.join(protein)


def find_orfs(seq, min_aa_length=50):
    """Basit ORF bulucu (doğrudan 3 okuma çerçevesi, + strand). Return list of dicts."""
    seq = seq.upper()
    orfs = []
    L = len(seq)
    for frame in range(3):
        i = frame
        while i < L - 2:
            codon = seq[i:i+3]
            if codon == 'ATG':  # start
                j = i
                while j < L - 2:
                    codon2 = seq[j:j+3]
                    if codon2 in ('TAA', 'TAG', 'TGA'):
                        prot_len = (j + 3 - i)//3
                        if prot_len >= min_aa_length:
                            dna_seq = seq[i:j+3]
                            prot = translate_dna(dna_seq)
                            orfs.append({'start': i+1, 'end': j+3, 'frame': frame+1, 'length_aa': prot_len, 'protein': prot})
                        i = j + 3
                        break
                    j += 3
                else:
                    # no stop codon found
                    break
            else:
                i += 3
    # Sort by length desc
    orfs.sort(key=lambda x: x['length_aa'], reverse=True)
    return orfs


# ------------------ Restriction Enzymes definitions ------------------
# Add more enzymes here if you like.
RESTRICTION_ENZYMES = {
    'EcoRI': 'GAATTC',
    'HindIII': 'AAGCTT',
    'BamHI': 'GGATCC',
    'NotI': 'GCGGCCGC',
    'XhoI': 'CTCGAG',
    'SmaI': 'CCCGGG',
    'PstI': 'CTGCAG',
    'NheI': 'GCTAGC',
    'KpnI': 'GGTACC',
    'SacI': 'GAGCTC',
    'SalI': 'GTCGAC'
}


# ---------------------- GUI ----------------------

class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=3, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super().__init__(fig)
        self.setParent(parent)

    def plot_bar(self, freq_dict):
        self.axes.clear()
        labels = ['A', 'C', 'G', 'T']
        values = [freq_dict.get(l, 0) for l in labels]
        self.axes.bar(labels, values)
        self.axes.set_title('Nucleotide counts')
        self.draw()

    def plot_gc(self, positions, gcs):
        self.axes.clear()
        if len(positions) == 0:
            self.axes.text(0.5, 0.5, 'No GC data', ha='center')
        else:
            self.axes.plot(positions, gcs)
            self.axes.set_title('Sliding-window GC%')
            self.axes.set_xlabel('Position')
            self.axes.set_ylabel('GC%')
        self.draw()


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Genetik PyQt5 Analiz Aracı - Prototip')
        self.resize(1000, 600)

        self.sequences = {}  # id -> seq

        self._create_menu()
        self._create_main_ui()

    def _create_menu(self):
        menubar = self.menuBar()
        file_menu = menubar.addMenu('File')

        open_action = QAction('Open FASTA', self)
        open_action.triggered.connect(self.load_fasta)
        file_menu.addAction(open_action)

        save_action = QAction('Export motif results (CSV)', self)
        save_action.triggered.connect(self.export_motif_results)
        file_menu.addAction(save_action)

        exit_action = QAction('Exit', self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

    def _create_main_ui(self):
        splitter = QSplitter(Qt.Horizontal)

        # Left: sequence list
        self.seq_list = QListWidget()
        self.seq_list.itemClicked.connect(self.on_seq_selected)
        splitter.addWidget(self.seq_list)

        # Right: tabs
        self.tabs = QTabWidget()

        # Sequence view tab
        seq_tab = QWidget()
        seq_layout = QVBoxLayout()
        self.seq_view = QTextEdit()
        self.seq_view.setReadOnly(True)
        seq_layout.addWidget(self.seq_view)
        seq_tab.setLayout(seq_layout)
        self.tabs.addTab(seq_tab, 'Sequence')

        # Analysis tab
        analysis_tab = QWidget()
        a_layout = QVBoxLayout()
        self.stats_label = QLabel('No sequence selected')
        a_layout.addWidget(self.stats_label)
        self.plot_canvas = PlotCanvas(self, width=5, height=3)
        a_layout.addWidget(self.plot_canvas)
        # sliding window controls
        sw_layout = QHBoxLayout()
        sw_layout.addWidget(QLabel('Sliding GC window:'))
        self.win_spin = QSpinBox()
        self.win_spin.setValue(50)
        self.win_spin.setMinimum(1)
        sw_layout.addWidget(self.win_spin)
        self.plot_gc_btn = QPushButton('Plot GC')
        self.plot_gc_btn.clicked.connect(self.on_plot_gc)
        sw_layout.addWidget(self.plot_gc_btn)
        a_layout.addLayout(sw_layout)
        analysis_tab.setLayout(a_layout)
        self.tabs.addTab(analysis_tab, 'Analysis')

        # Motif tab
        motif_tab = QWidget()
        m_layout = QVBoxLayout()
        input_layout = QHBoxLayout()
        input_layout.addWidget(QLabel('Motif (or regex):'))
        self.motif_input = QLineEdit()
        input_layout.addWidget(self.motif_input)
        self.regex_checkbox = QCheckBox('Regex')
        input_layout.addWidget(self.regex_checkbox)
        self.revcomp_checkbox = QCheckBox('Search reverse-complement')
        input_layout.addWidget(self.revcomp_checkbox)
        self.search_btn = QPushButton('Search motif in all sequences')
        self.search_btn.clicked.connect(self.on_search_motif)
        input_layout.addWidget(self.search_btn)
        m_layout.addLayout(input_layout)

        self.motif_table = QTableWidget(0, 6)
        self.motif_table.setHorizontalHeaderLabels(['SeqID', 'Start', 'End', 'Strand', 'Match', 'Context'])
        self.motif_table.cellClicked.connect(self.on_motif_cell_clicked)
        m_layout.addWidget(self.motif_table)

        motif_tab.setLayout(m_layout)
        self.tabs.addTab(motif_tab, 'Motifs')

        # ORF tab
        orf_tab = QWidget()
        orf_layout = QVBoxLayout()
        orf_ctrl_layout = QHBoxLayout()
        orf_ctrl_layout.addWidget(QLabel('Minimum AA length:'))
        self.min_orf_spin = QSpinBox()
        self.min_orf_spin.setMinimum(10)
        self.min_orf_spin.setValue(50)
        orf_ctrl_layout.addWidget(self.min_orf_spin)
        self.find_orf_btn = QPushButton('Find ORFs (current sequence)')
        self.find_orf_btn.clicked.connect(self.on_find_orfs)
        orf_ctrl_layout.addWidget(self.find_orf_btn)
        orf_layout.addLayout(orf_ctrl_layout)
        self.orf_table = QTableWidget(0, 5)
        self.orf_table.setHorizontalHeaderLabels(['Start', 'End', 'Frame', 'AA Length', 'Protein (first 50 aa)'])
        self.orf_table.cellClicked.connect(self.on_orf_cell_clicked)
        orf_layout.addWidget(self.orf_table)
        orf_tab.setLayout(orf_layout)
        self.tabs.addTab(orf_tab, 'ORFs')

        # ------------------ Restriction Enzymes tab ------------------
        restr_tab = QWidget()
        restr_layout = QVBoxLayout()

        # controls: enzyme select, search all checkbox, button
        ctrl_layout = QHBoxLayout()
        ctrl_layout.addWidget(QLabel('Enzyme:'))
        self.enzyme_combo = QComboBox()
        # populate combo with enzyme names
        for name in sorted(RESTRICTION_ENZYMES.keys()):
            self.enzyme_combo.addItem(name)
        ctrl_layout.addWidget(self.enzyme_combo)

        self.search_all_checkbox = QCheckBox('Search all sequences (if unchecked, search current)')
        ctrl_layout.addWidget(self.search_all_checkbox)

        self.find_restr_btn = QPushButton('Find restriction sites')
        self.find_restr_btn.clicked.connect(self.on_find_restriction_sites)
        ctrl_layout.addWidget(self.find_restr_btn)

        restr_layout.addLayout(ctrl_layout)

        # results table
        self.restr_table = QTableWidget(0, 5)
        self.restr_table.setHorizontalHeaderLabels(['SeqID', 'Enzyme', 'Position', 'Strand', 'Site'])
        self.restr_table.cellClicked.connect(self.on_restr_cell_clicked)
        restr_layout.addWidget(self.restr_table)

        restr_tab.setLayout(restr_layout)
        self.tabs.addTab(restr_tab, 'Restriction Enzymes')
        # -----------------------------------------------------------

        splitter.addWidget(self.tabs)
        splitter.setStretchFactor(1, 3)

        self.setCentralWidget(splitter)

        # motif results cache for export
        self._motif_results = []

    # ------------------ UI callbacks ------------------

    def load_fasta(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Open FASTA', '', 'FASTA files (*.fasta *.fa *.fna);;All files (*)')
        if not path:
            return
        try:
            seqs = parse_fasta(path)
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'FASTA parse error: {e}')
            return
        self.sequences = seqs
        self.seq_list.clear()
        for sid in self.sequences.keys():
            self.seq_list.addItem(sid)
        self.statusBar().showMessage(f'Loaded {len(self.sequences)} sequences from {os.path.basename(path)}')

    def on_seq_selected(self, item):
        sid = item.text()
        seq = self.sequences.get(sid, '')
        self.display_sequence(seq)
        self.update_stats_display(sid, seq)

    def display_sequence(self, seq):
        # Show sequence with monospace-like formatting
        self.seq_view.setPlainText(seq)

    def update_stats_display(self, sid, seq):
        freqs, total = nucleotide_freq(seq)
        gc = gc_content(seq)
        stats = f"SeqID: {sid}\nLength: {len(seq)} nt\nGC%: {gc:.2f}\nA: {freqs['A']} C: {freqs['C']} G: {freqs['G']} T: {freqs['T']}"
        self.stats_label.setText(stats)
        self.plot_canvas.plot_bar(freqs)

    def on_plot_gc(self):
        cur_item = self.seq_list.currentItem()
        if not cur_item:
            QMessageBox.information(self, 'Info', 'Önce bir sekans seçiniz.')
            return
        sid = cur_item.text()
        seq = self.sequences.get(sid, '')
        w = int(self.win_spin.value())
        pos, gcs = sliding_gc(seq, window=w)
        self.plot_canvas.plot_gc(pos, gcs)

    def on_search_motif(self):
        motif = self.motif_input.text().strip()
        if not motif:
            QMessageBox.information(self, 'Info', 'Motif girin.')
            return
        use_regex = self.regex_checkbox.isChecked()
        search_rev = self.revcomp_checkbox.isChecked()
        self.motif_table.setRowCount(0)
        results = []
        for sid, seq in self.sequences.items():
            s = seq.upper()
            # forward
            if use_regex:
                try:
                    for m in re.finditer(motif, s):
                        start = m.start() + 1
                        end = m.end()
                        match = m.group()
                        context = s[max(0, start-6):min(len(s), end+5)]
                        results.append((sid, start, end, '+', match, context))
                except re.error as e:
                    QMessageBox.critical(self, 'Regex error', f'Regex hatası: {e}')
                    return
            else:
                # simple substring search
                sub = motif.upper()
                idx = s.find(sub)
                while idx != -1:
                    start = idx + 1
                    end = idx + len(sub)
                    match = s[idx:idx+len(sub)]
                    context = s[max(0, idx-5):min(len(s), idx+len(sub)+5)]
                    results.append((sid, start, end, '+', match, context))
                    idx = s.find(sub, idx+1)

            # reverse complement
            if search_rev:
                rc = revcomp(s)
                if use_regex:
                    try:
                        for m in re.finditer(motif, rc):
                            # convert positions back to original
                            rc_start = m.start()
                            rc_end = m.end() - 1
                            # mapping: pos_in_orig = len(s) - pos_in_rc
                            # start_orig = len(s) - rc_end
                            start = len(s) - (m.end() - 1)
                            end = len(s) - m.start()
                            match = m.group()
                            context = s[max(0, start-6):min(len(s), end+5)]
                            results.append((sid, start, end, '-', match, context))
                    except re.error as e:
                        QMessageBox.critical(self, 'Regex error', f'Regex hatası: {e}')
                        return
                else:
                    sub = motif.upper()
                    idx = rc.find(sub)
                    while idx != -1:
                        # convert to original coordinates
                        start = len(s) - (idx + len(sub)) + 1
                        end = len(s) - idx
                        match = s[start-1:end]
                        context = s[max(0, start-6):min(len(s), end+5)]
                        results.append((sid, start, end, '-', match, context))
                        idx = rc.find(sub, idx+1)

        # populate table
        self._motif_results = []
        for row, r in enumerate(results):
            self.motif_table.insertRow(row)
            for col, val in enumerate(r):
                self.motif_table.setItem(row, col, QTableWidgetItem(str(val)))
            self._motif_results.append({'seqid': r[0], 'start': r[1], 'end': r[2], 'strand': r[3], 'match': r[4]})
        self.statusBar().showMessage(f'Found {len(results)} motif matches')

    def on_motif_cell_clicked(self, row, col):
        # highlight match in sequence viewer if the seq is currently displayed
        sid_item = self.motif_table.item(row, 0)
        if sid_item is None:
            return
        sid = sid_item.text()
        cur_item = self.seq_list.currentItem()
        if cur_item and cur_item.text() != sid:
            # select the row in seq_list
            items = self.seq_list.findItems(sid, Qt.MatchExactly)
            if items:
                self.seq_list.setCurrentItem(items[0])
                self.display_sequence(self.sequences[sid])
                self.update_stats_display(sid, self.sequences[sid])
        start = int(self.motif_table.item(row, 1).text())
        end = int(self.motif_table.item(row, 2).text())
        self.highlight_region(start-1, end-1)

    def highlight_region(self, start_idx, end_idx):
        # clear formatting
        cursor = self.seq_view.textCursor()
        fmt_clear = QTextCharFormat()
        fmt_clear.setBackground(QColor('white'))
        cursor.select(QTextCursor.Document)
        cursor.setCharFormat(fmt_clear)

        # highlight desired region
        cursor = self.seq_view.textCursor()
        cursor.setPosition(0)
        cursor.movePosition(QTextCursor.Right, QTextCursor.MoveAnchor, start_idx)
        cursor.movePosition(QTextCursor.Right, QTextCursor.KeepAnchor, end_idx - start_idx + 1)
        fmt = QTextCharFormat()
        fmt.setBackground(QColor('yellow'))
        cursor.setCharFormat(fmt)
        self.seq_view.setTextCursor(cursor)

    def export_motif_results(self):
        if not self._motif_results:
            QMessageBox.information(self, 'Info', 'Export edilecek motif sonucu yok.')
            return
        path, _ = QFileDialog.getSaveFileName(self, 'Save motif results', 'motif_results.csv', 'CSV files (*.csv)')
        if not path:
            return
        try:
            with open(path, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['seqid', 'start', 'end', 'strand', 'match'])
                for r in self._motif_results:
                    writer.writerow([r['seqid'], r['start'], r['end'], r['strand'], r['match']])
            QMessageBox.information(self, 'Saved', f'Saved {len(self._motif_results)} rows to {os.path.basename(path)}')
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'CSV yazma hatası: {e}')

    def on_find_orfs(self):
        cur_item = self.seq_list.currentItem()
        if not cur_item:
            QMessageBox.information(self, 'Info', 'Önce bir sekans seçiniz.')
            return
        sid = cur_item.text()
        seq = self.sequences.get(sid, '')
        min_len = int(self.min_orf_spin.value())
        orfs = find_orfs(seq, min_aa_length=min_len)
        self.orf_table.setRowCount(0)
        for i, o in enumerate(orfs):
            self.orf_table.insertRow(i)
            self.orf_table.setItem(i, 0, QTableWidgetItem(str(o['start'])))
            self.orf_table.setItem(i, 1, QTableWidgetItem(str(o['end'])))
            self.orf_table.setItem(i, 2, QTableWidgetItem(str(o['frame'])))
            self.orf_table.setItem(i, 3, QTableWidgetItem(str(o['length_aa'])))
            prot_preview = o['protein'][:50]
            self.orf_table.setItem(i, 4, QTableWidgetItem(prot_preview))
        self.statusBar().showMessage(f'Found {len(orfs)} ORFs in {sid}')
        # cache last ORFs
        self._last_orfs = orfs

    def on_orf_cell_clicked(self, row, col):
        if not hasattr(self, '_last_orfs'):
            return
        o = self._last_orfs[row]
        # show full protein in a message box or new window
        msg = QMessageBox(self)
        msg.setWindowTitle(f"Protein: aa {o['length_aa']} (frame {o['frame']})")
        prot_text = o['protein']
        # limit size in message box
        msg.setText(prot_text[:5000])
        msg.exec_()

    # ---------------- Restriction enzyme search -----------------

    def on_find_restriction_sites(self):
        enzyme_name = self.enzyme_combo.currentText()
        if not enzyme_name:
            QMessageBox.information(self, 'Info', 'Enzim seçin.')
            return
        pattern = RESTRICTION_ENZYMES.get(enzyme_name, '').upper()
        if not pattern:
            QMessageBox.critical(self, 'Error', 'Seçilen enzimin dizisi bulunamadı.')
            return

        search_all = self.search_all_checkbox.isChecked()
        self.restr_table.setRowCount(0)
        row = 0

        targets = []
        if search_all:
            targets = list(self.sequences.items())
        else:
            cur_item = self.seq_list.currentItem()
            if not cur_item:
                QMessageBox.information(self, 'Info', 'Önce bir sekans seçiniz (veya "Search all" işaretleyin).')
                return
            sid = cur_item.text()
            targets = [(sid, self.sequences.get(sid, ''))]

        # search forward pattern and reverse complement of pattern
        pattern_rc = revcomp(pattern)

        for sid, seq in targets:
            s = seq.upper()
            # forward occurrences
            start_idx = 0
            while True:
                idx = s.find(pattern, start_idx)
                if idx == -1:
                    break
                self.restr_table.insertRow(row)
                self.restr_table.setItem(row, 0, QTableWidgetItem(sid))
                self.restr_table.setItem(row, 1, QTableWidgetItem(enzyme_name))
                self.restr_table.setItem(row, 2, QTableWidgetItem(str(idx+1)))
                self.restr_table.setItem(row, 3, QTableWidgetItem('+'))
                self.restr_table.setItem(row, 4, QTableWidgetItem(pattern))
                row += 1
                start_idx = idx + 1  # allow overlapping

            # reverse complement occurrences (pattern on - strand)
            start_idx = 0
            # If pattern == pattern_rc and we already found forward, we will also find here duplicates;
            # that's fine, but we can avoid duplicates by checking equality:
            if pattern_rc != pattern:
                while True:
                    idx = s.find(pattern_rc, start_idx)
                    if idx == -1:
                        break
                    self.restr_table.insertRow(row)
                    self.restr_table.setItem(row, 0, QTableWidgetItem(sid))
                    self.restr_table.setItem(row, 1, QTableWidgetItem(enzyme_name))
                    self.restr_table.setItem(row, 2, QTableWidgetItem(str(idx+1)))
                    self.restr_table.setItem(row, 3, QTableWidgetItem('-'))
                    self.restr_table.setItem(row, 4, QTableWidgetItem(pattern_rc))
                    row += 1
                    start_idx = idx + 1

        self.statusBar().showMessage(f'Restriction sites found: {row}')

    def on_restr_cell_clicked(self, row, col):
        # highlight site in sequence viewer when clicked
        sid_item = self.restr_table.item(row, 0)
        pos_item = self.restr_table.item(row, 2)
        site_item = self.restr_table.item(row, 4)
        if not sid_item or not pos_item or not site_item:
            return
        sid = sid_item.text()
        pos = int(pos_item.text())
        site = site_item.text()
        # select sequence in left list if different
        cur_item = self.seq_list.currentItem()
        if not cur_item or cur_item.text() != sid:
            items = self.seq_list.findItems(sid, Qt.MatchExactly)
            if items:
                self.seq_list.setCurrentItem(items[0])
                self.display_sequence(self.sequences[sid])
                self.update_stats_display(sid, self.sequences[sid])
        start_idx = pos - 1
        end_idx = start_idx + len(site) - 1
        self.highlight_region(start_idx, end_idx)

if __name__ == '__main__':
    app = QApplication(sys.argv)

    import qdarkstyle

    app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())

    win = MainWindow()
    win.show()
    sys.exit(app.exec_())

