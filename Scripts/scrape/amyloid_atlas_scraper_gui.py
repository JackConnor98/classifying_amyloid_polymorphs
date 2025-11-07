import requests
from bs4 import BeautifulSoup
import pandas as pd
import tkinter as tk
from tkinter import ttk, messagebox
import os

URL = "https://people.mbi.ucla.edu/sawaya/amyloidatlas/"

# -------------------------------------------------------------------
# DATA FETCHING
# -------------------------------------------------------------------
def fetch_html(url):
    response = requests.get(url)
    return BeautifulSoup(response.text, "html.parser")

def parse_table_data(doc):
    table_values = doc.find_all("td")
    headers = []
    rows = []
    current_row = []

    for i in range(len(table_values)):
        if i < 9 and i % 9 not in (5, 6):
            headers.append(table_values[i].get_text(strip=True))

        if i >= 9 and i % 9 not in (5, 6):
            current_row.append(table_values[i].get_text(strip=True))

        if i > 8 and i % 9 == 8:
            rows.append(current_row)
            current_row = []

    return headers, rows

# -------------------------------------------------------------------
# GUI APPLICATION
# -------------------------------------------------------------------
class App(tk.Tk):
    def __init__(self, df):
        super().__init__()
        self.title("Amyloid Atlas Scraper")
        self.geometry("1200x700")
        self.df = df
        self.selected_rows = set()  # track highlighted rows
        self.create_widgets()
        self.populate_table(self.df)

    # ---------------------------------------------------------------
    def create_widgets(self):
        # Entry + buttons frame
        entry_frame = tk.Frame(self)
        entry_frame.pack(pady=10)

        tk.Label(entry_frame, text="Enter Protein Pattern:").pack(side="left", padx=5)
        self.pattern_entry = tk.Entry(entry_frame, width=30)
        self.pattern_entry.pack(side="left", padx=5)

        ttk.Button(entry_frame, text="Filter", command=self.filter_table).pack(side="left", padx=5)
        ttk.Button(entry_frame, text="Reset", command=self.reset_table).pack(side="left", padx=5)
        ttk.Button(entry_frame, text="Select All", command=self.select_all).pack(side="left", padx=5)
        ttk.Button(entry_frame, text="Clear All", command=self.clear_all).pack(side="left", padx=5)
        ttk.Button(entry_frame, text="Export Selected", command=self.export_selected).pack(side="left", padx=5)

        # Table frame
        table_frame = tk.Frame(self)
        table_frame.pack(fill="both", expand=True)

        # Scrollbars
        yscroll = tk.Scrollbar(table_frame, orient="vertical")
        yscroll.pack(side="right", fill="y")
        xscroll = tk.Scrollbar(table_frame, orient="horizontal")
        xscroll.pack(side="bottom", fill="x")

        # Treeview (table)
        self.tree = ttk.Treeview(
            table_frame,
            columns=list(self.df.columns),
            show="headings",
            yscrollcommand=yscroll.set,
            xscrollcommand=xscroll.set,
            selectmode="none"
        )

        yscroll.config(command=self.tree.yview)
        xscroll.config(command=self.tree.xview)

        for col in self.df.columns:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=150, anchor="center")

        self.tree.pack(fill="both", expand=True)
        self.tree.bind("<Button-1>", self.on_row_click)

        # Tag style for highlighted rows
        self.tree.tag_configure("selected", background="#cce5ff")  # light blue

    # ---------------------------------------------------------------
    def populate_table(self, df):
        self.tree.delete(*self.tree.get_children())
        self.selected_rows.clear()
        for _, row in df.iterrows():
            row_id = self.tree.insert("", "end", values=list(row))

    # ---------------------------------------------------------------
    def filter_table(self):
        pattern = self.pattern_entry.get().strip()
        if not pattern:
            messagebox.showwarning("Warning", "Enter a protein pattern first.")
            return
        filtered = self.df[self.df["Protein"].str.contains(pattern, case=False, na=False)]
        if filtered.empty:
            messagebox.showinfo("No Results", f"No entries match pattern '{pattern}'.")
            return
        self.populate_table(filtered)

    # ---------------------------------------------------------------
    def reset_table(self):
        self.pattern_entry.delete(0, tk.END)
        self.populate_table(self.df)

    # ---------------------------------------------------------------
    def clear_all(self):
        for row_id in self.tree.get_children():
            self.tree.item(row_id, tags=())
        self.selected_rows.clear()

    def select_all(self):
        for row_id in self.tree.get_children():
            self.tree.item(row_id, tags=("selected",))
            self.selected_rows.add(row_id)

    # ---------------------------------------------------------------
    def on_row_click(self, event):
        row_id = self.tree.identify_row(event.y)
        if not row_id:
            return
        # toggle highlight
        if row_id in self.selected_rows:
            self.selected_rows.remove(row_id)
            self.tree.item(row_id, tags=())
        else:
            self.selected_rows.add(row_id)
            self.tree.item(row_id, tags=("selected",))

    # ---------------------------------------------------------------
    def export_selected(self):
        if not self.selected_rows:
            messagebox.showwarning("Nothing Selected", "Select something first, noob.")
            return

        # Gather selected rows
        selected_data = [self.tree.item(row_id, "values") for row_id in self.selected_rows]
        export_df = pd.DataFrame(selected_data, columns=self.df.columns)

        # Filter out non-CryoEM rows before saving
        cryoEM_df = export_df[export_df["Method"].str.lower() == "cryoem"]
        non_cryoEM_df = export_df[export_df["Method"].str.lower() != "cryoem"]

        # --- Save excluded non-CryoEM entries ---
        excluded_file = os.path.join("Output", "excluded_list.txt")
        os.makedirs("Output", exist_ok=True)
        file_exists = os.path.exists(excluded_file)
        existing_pdbs = set()
        if file_exists:
            with open(excluded_file, "r") as f:
                for line in f:
                    if line.strip() and not line.startswith("PDB"):
                        existing_pdbs.add(line.split("\t")[0])

        with open(excluded_file, "a") as f:
            if not file_exists:
                f.write("PDB\tReason\n")
            for _, row in non_cryoEM_df.iterrows():
                pdb = row["PDB ID"]
                if pdb not in existing_pdbs:
                    method = row["Method"]
                    f.write(f"{pdb}\tNot solved using CryoEM, solved using {method}\n")

        # --- Save only CryoEM entries ---
        cryoEM_df.to_csv(os.path.join("Output", "selected_pdbs_metadata.csv"), index=False)

        # Save only PDB IDs
        cryoEM_df.rename(columns={"PDB ID": "PDB"})["PDB"].to_csv(
            os.path.join("Output", "pdb_names.txt"), index=False, header=True
        )

        messagebox.showinfo(
            "Success",
            f"Saved {len(cryoEM_df)} entries; {len(non_cryoEM_df)} excluded"
        )


# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------
def main():
    print("Fetching data...")
    doc = fetch_html(URL)
    headers, rows = parse_table_data(doc)
    df = pd.DataFrame(rows, columns=headers)

    # Fix Resolution column name
    for col in df.columns:
        if "resol" in col.lower():
            df.rename(columns={col: "Resolution"}, inplace=True)
            break

    app = App(df)
    app.mainloop()


if __name__ == "__main__":
    main()
