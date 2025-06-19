import sqlite3
import tkinter as tk
from tkinter import ttk
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def fetch_data(db_path, table_name, order_by=None):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    if order_by:
        cursor.execute(f"SELECT * FROM {table_name} ORDER BY {order_by} ASC")
    else:
        cursor.execute(f"SELECT * FROM {table_name}")
    columns = [description[0] for description in cursor.description]
    rows = cursor.fetchall()
    conn.close()
    return columns, rows

def parse_xyz(xyz_text):
    lines = xyz_text.strip().splitlines()
    atoms = []
    coords = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) == 4:
            atoms.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return atoms, coords

def draw_xyz_3d_on_canvas(xyz_text, canvas_frame, title="3D-Struktur"):
    import numpy as np
    for widget in canvas_frame.winfo_children():
        widget.destroy()
    atoms, coords = parse_xyz(xyz_text)
    xs, ys, zs = zip(*coords)
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xs, ys, zs, s=100, c='b', alpha=0.7)
    for i, atom in enumerate(atoms):
        ax.text(xs[i], ys[i], zs[i], atom, size=10, zorder=1, color='k')
    # Linien zu jeweils nächstem Nachbarn zeichnen
    coords_np = np.array(coords)
    drawn = set()
    for i, c1 in enumerate(coords_np):
        dists = np.linalg.norm(coords_np - c1, axis=1)
        dists[i] = np.inf  # sich selbst ignorieren
        j = np.argmin(dists)
        key = tuple(sorted((i, j)))
        if key not in drawn:
            ax.plot([c1[0], coords_np[j][0]], [c1[1], coords_np[j][1]], [c1[2], coords_np[j][2]], c='gray', lw=2)
            drawn.add(key)
    ax.set_title(title)
    # Koordinatensystem wieder anzeigen (keine set_axis_off, keine Ticks ausblenden)
    canvas = FigureCanvasTkAgg(fig, master=canvas_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(expand=True, fill='both')
    plt.close(fig)

def show_data():
    root = tk.Tk()
    root.title("Datenbank-Viewer")

    notebook = ttk.Notebook(root)
    notebook.pack(expand=True, fill='both')

    # Metalle-Tab
    metalle_frame = ttk.Frame(notebook)
    notebook.add(metalle_frame, text="Metals")

    metalle_db_path = os.path.join(os.path.dirname(__file__), "metals.db")
    metalle_columns, metalle_rows = fetch_data(metalle_db_path, "metalle", order_by="ordnungszahl")

    # German to English mapping for metals columns
    metals_column_translation = {
        "id": "ID",
        "name": "Name",
        "ordnungszahl": "Atomic Number",
        "d_elektronen": "d Electrons",
        "oxidation": "Oxidation State",
        "koordinationszahl": "Coordination Number",
        "geometrie": "Geometry"
    }

    metalle_tree = ttk.Treeview(metalle_frame, columns=metalle_columns, show='headings')
    for col in metalle_columns:
        metalle_tree.heading(col, text=metals_column_translation.get(col, col))
        metalle_tree.column(col, width=120)
    for row in metalle_rows:
        metalle_tree.insert('', tk.END, values=row)
    metalle_tree.pack(expand=True, fill='both')

    # Liganden-Tab
    liganden_frame = ttk.Frame(notebook)
    notebook.add(liganden_frame, text="Ligands")

    liganden_db_path = os.path.join(os.path.dirname(__file__), "ligands.db")
    liganden_columns, liganden_rows = fetch_data(liganden_db_path, "liganden")

    # German to English mapping for ligands columns
    ligands_column_translation = {
        "id": "ID",
        "name": "Name",
        "ladung": "Charge",
        "xyz_daten": "XYZ Data"
    }

    display_columns = liganden_columns[:-1]
    liganden_tree = ttk.Treeview(liganden_frame, columns=display_columns, show='headings')
    for col in display_columns:
        liganden_tree.heading(col, text=ligands_column_translation.get(col, col))
        liganden_tree.column(col, width=120)
    for row in liganden_rows:
        liganden_tree.insert('', tk.END, values=row[:-1])
    liganden_tree.pack(side='left', expand=True, fill='both')

    # Frame für 3D-Plot
    plot_frame = tk.Frame(liganden_frame)
    plot_frame.pack(side='right', expand=True, fill='both')

    def on_ligand_select(event=None):
        selected = liganden_tree.selection()
        if not selected and liganden_rows:
            idx = 0
        elif selected:
            idx = liganden_tree.index(selected[0])
        else:
            return
        xyz_text = liganden_rows[idx][-1]
        name = liganden_rows[idx][1]
        draw_xyz_3d_on_canvas(xyz_text, plot_frame, title=name)

    liganden_tree.bind('<<TreeviewSelect>>', on_ligand_select)

    # Direkt beim Start die erste Struktur anzeigen
    if liganden_rows:
        liganden_tree.selection_set(liganden_tree.get_children()[0])
        on_ligand_select()

    root.mainloop()

if __name__ == "__main__":
    show_data()