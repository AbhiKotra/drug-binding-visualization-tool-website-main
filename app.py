from flask import Flask, request, jsonify, render_template
import os
import subprocess
import base64
import uuid
import sys
import csv
from datetime import datetime

# ─────────────────────────────────────────────────────────────────
# app.py — the Flask web server
#
# Think of this file as the "manager" of your restaurant.
# It receives files from the browser, tells the kitchen (your scripts)
# to run, then sends results back to the browser.
# ─────────────────────────────────────────────────────────────────

app = Flask(__name__)

# Where this file lives (project root)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Folders for temporary uploads and results
UPLOAD_FOLDER  = os.path.join(BASE_DIR, "data", "uploads")
RESULTS_FOLDER = os.path.join(BASE_DIR, "data", "results")
IMAGES_FOLDER  = os.path.join(BASE_DIR, "data", "images")
SCRIPTS_FOLDER = os.path.join(BASE_DIR, "scripts")

# Create these folders if they don't exist yet
os.makedirs(UPLOAD_FOLDER,  exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)
os.makedirs(IMAGES_FOLDER,  exist_ok=True)


# ── Route 1: Serve the main webpage ───────────────────────────────
# When someone visits your URL, this sends them index.html
@app.route('/')
def index():
    return render_template('index.html')


# ── Route 2: Receive uploaded files & run analysis ────────────────
# When the user clicks Analyze, the browser sends the files here
@app.route('/analyze', methods=['POST'])
def analyze():
    # ── Step 1: Get the uploaded files ──
    if 'alphafold' not in request.files or 'bound' not in request.files:
        return jsonify({'error': 'Please upload both PDB files.'}), 400

    alphafold_file = request.files['alphafold']
    bound_file     = request.files['bound']
    ligand_code    = request.form.get('ligand_code', 'IRE').strip().upper()
    protein_name = request.form.get('protein_name', 'Unknown Protein').strip()
    mutation_name = request.form.get('mutation_name', 'Unknown Mutation').strip()
    disease_name = request.form.get('disease_name', 'Unknown Disease').strip()

    # ── Step 2: Save with unique names so simultaneous users don't overwrite each other ──
    session_id     = str(uuid.uuid4())[:8]   # random 8-character ID
    alphafold_path = os.path.join(UPLOAD_FOLDER, f"{session_id}_alphafold.pdb")
    bound_path     = os.path.join(UPLOAD_FOLDER, f"{session_id}_bound.pdb")
    results_dir    = os.path.join(RESULTS_FOLDER, session_id)
    images_dir     = os.path.join(IMAGES_FOLDER,  session_id)

    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(images_dir,  exist_ok=True)

    alphafold_file.save(alphafold_path)
    bound_file.save(bound_path)
    # Save analysis to history
    history_path = os.path.join(BASE_DIR, "data", "history.csv")
    history_exists = os.path.exists(history_path)

    with open(history_path, "a", newline="") as f:
        writer = csv.writer(f)

        if not history_exists:
            writer.writerow([
                "session_id",
                "protein_name",
                "mutation_name",
                "disease_name",
                "ligand_code",
                "alphafold_file",
                "bound_file",
                "timestamp"
            ])

        writer.writerow([
            session_id,
            protein_name,
            mutation_name,
            disease_name,
            ligand_code,
            alphafold_file.filename,
            bound_file.filename,
            datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        ])

    python = sys.executable  # use the same Python that's running Flask

    # ── Step 3: Run analyze_structure.py ──
    try:
        subprocess.run([
            python,
            os.path.join(SCRIPTS_FOLDER, 'analyze_structure.py'),
            alphafold_path,
            bound_path,
            results_dir,
            ligand_code
        ], check=True, capture_output=True, text=True)

        subprocess.run([
            python,
            os.path.join(SCRIPTS_FOLDER, 'find_binding_pocket.py'),
            bound_path,
            ligand_code,
            results_dir
        ], check=True, capture_output=True, text=True)

        subprocess.run([
            python,
            os.path.join(SCRIPTS_FOLDER, 'plot_summary.py'),
            results_dir,
            images_dir
        ], check=True, capture_output=True, text=True)

    except subprocess.CalledProcessError as e:
        return jsonify({
            'error': 'Backend script failed',
            'details': e.stderr
        }), 500
    # ── Step 6: Read the results and send them back to the browser ──

    # Read the CSV stats
    csv_path = os.path.join(results_dir, "structure_summary.csv")
    stats = []
    if os.path.exists(csv_path):
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            stats = list(reader)
            
    advanced_path = os.path.join(results_dir, "advanced_metrics.csv")
    advanced_metrics = []
    if os.path.exists(advanced_path):
        with open(advanced_path) as f:
            reader = csv.DictReader(f)
            advanced_metrics = list(reader)

    # Read the binding pocket residues text file
    pocket_path = os.path.join(results_dir, "binding_pocket_residues.txt")
    pocket_text = ""
    if os.path.exists(pocket_path):
        with open(pocket_path) as f:
            pocket_text = f.read()

    # Convert the graph PNG to base64 (a way to embed an image directly in JSON)
    graphs = []
    for filename in os.listdir(images_dir):
        if filename.endswith('.png'):
            filepath = os.path.join(images_dir, filename)
            with open(filepath, 'rb') as f:
                encoded = base64.b64encode(f.read()).decode('utf-8')
                graphs.append({'name': filename, 'data': encoded})

    # Send everything back to the browser as JSON
    return jsonify({
        'status': 'success',
        'stats': stats,
        'advanced_metrics': advanced_metrics,
        'pocket': pocket_text,
        'graphs': graphs,
        'alphafold_id': f"{session_id}_alphafold",
        'bound_id': f"{session_id}_bound",
        'results_folder': session_id,
        'metadata': {
            'protein_name': protein_name,
            'mutation_name': mutation_name,
            'ligand_code': ligand_code,
            'disease_name': disease_name
        },
})


# ── Route 3: Serve uploaded PDB files to the 3D viewer ────────────
# The NGL viewer in the browser needs to fetch the PDB file as a URL
@app.route('/pdb/<filename>')
def serve_pdb(filename):
    from flask import send_from_directory
    return send_from_directory(UPLOAD_FOLDER, filename + '.pdb',
                               mimetype='text/plain')

@app.route('/download/<session_id>/<filename>')
def download_file(session_id, filename):

    from flask import send_from_directory

    folder = os.path.join(RESULTS_FOLDER, session_id)

    return send_from_directory(folder, filename, as_attachment=True)
# ── Start the server ──────────────────────────────────────────────
@app.route('/history')
def history():
    history_path = os.path.join(BASE_DIR, "data", "history.csv")

    rows = []

    if os.path.exists(history_path):
        with open(history_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

    rows.reverse()

    return jsonify(rows)

if __name__ == '__main__':
    app.run(debug=False, host='0.0.0.0', port=int(os.environ.get('PORT', 5000)))
