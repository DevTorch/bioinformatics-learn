from pathlib import Path

from Bio import SeqIO, SearchIO
from Bio.Blast import NCBIWWW
from Bio.PDB import PDBIO, MMCIFIO, PDBParser, MMCIFParser

import nglview as nv
print(nv.__version__)

import warnings

import pymol
from pymol import cmd

from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.filterwarnings("ignore")
warnings.filterwarnings('ignore', category=PDBConstructionWarning)

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
fasta_file = DATA_DIR / 'longest_orf.fasta'
pdb_file = DATA_DIR / '6YYT.pdb'
output_html = DATA_DIR / 'structure_view.html'

# BLAST search
covid_protein_seq = SeqIO.read(fasta_file, 'fasta')
result_handle = NCBIWWW.qblast("blastp", "pdb", covid_protein_seq.seq)  # PDB – Protein Data Bank
blast_record = SearchIO.read(result_handle, "blast-xml")
print(blast_record[0:10])

for record in blast_record:
    print(f'Sequence ID: {record.id}')
    print(f'Description: {record.description}')
    print(f'E Value: {record[0].evalue}')
    print(f'Bit Score: {record[0].bitscore}')
    print(f'Alignment length: {record[0].aln}')

parser = PDBParser(QUIET=True)
structure = parser.get_structure('6YYT', pdb_file)
print(structure)

for chain in structure.get_chains():
    print(chain.get_id())

view = nv.show_biopython(structure)

pymol.finish_launching()

cmd.load("data/6YYT.pdb", "rdRp")
cmd.remove("solvent")
cmd.hide("everything")

# --- БАЗОВЫЕ СЕЛЕКТЫ ---
cmd.select("nsp12", "rdRp and chain A")
cmd.select("nsp8",  "rdRp and (chain B or chain D)")
cmd.select("nsp7",  "rdRp and chain C")
cmd.select("rna",   "rdRp and (chain P or chain Q or chain T or chain U)")

# --- ВИДЫ ---
cmd.show("cartoon", "nsp12 or nsp8 or nsp7")
cmd.show("sticks", "rna")
cmd.set("stick_radius", 0.22, "rna")

# --- ЦВЕТА (можешь поменять под вкус) ---
cmd.color("wheat",    "nsp12")
cmd.color("palecyan", "nsp7")
cmd.color("lightpink","nsp8")
cmd.color("gray70",   "rna")

# --- АККУРАТНОСТИ ВНЕШНЕГО ВИДА ---
cmd.bg_color("white")
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_smooth_loops", 1)
cmd.set("cartoon_transparency", 0.0)

# --- ПОДСВЕТКА АКТИВНОГО САЙТА (часто D760/D761 в nsp12) ---
cmd.select("active_site", "nsp12 and resi 760+761")
cmd.show("sticks", "active_site")
cmd.color("hotpink", "active_site")
cmd.label("active_site and name CA", '"%s%s" % (resn, resi)')
cmd.set("label_size", 18)

# --- ИНТЕРФЕЙС nsp12 ↔ RNA (в радиусе 4Å) ---
cmd.select("rna_iface", "byres (rna within 4 of nsp12)")
cmd.show("sticks", "rna_iface")
cmd.color("tv_orange", "rna_iface")

# --- ИНТЕРФЕЙСЫ nsp12 ↔ nsp8 (в радиусе 4Å) ---
cmd.select("iface_12", "byres (nsp12 within 4 of nsp8)")
cmd.select("iface_8",  "byres (nsp8 within 4 of nsp12)")
cmd.show("sticks", "iface_12 or iface_8")
cmd.color("orange", "iface_12")
cmd.color("marine", "iface_8")

cmd.orient("rdRp")

# --- КРАСИВЫЙ РЕНДЕР В PNG (опционально) ---
cmd.set("ray_opaque_background", 0)
cmd.set("antialias", 2)
cmd.ray(1800, 1400)
cmd.png("data/6YYT_preset.png", dpi=300)
print("Saved: data/6YYT_preset.png")