from flask import Flask, request, render_template, send_file, redirect, url_for
import subprocess
import sys
from pathlib import Path
from datetime import datetime
import uuid
import shutil

app = Flask(__name__, template_folder="templates", static_folder="static")

SCRIPT_PATH = Path(__file__).resolve().parent.parent / "comborate.py"
BASE_DIR = Path(__file__).resolve().parent
UPLOAD_DIR = BASE_DIR / "uploads"
UPLOAD_DIR.mkdir(exist_ok=True)


@app.route("/", methods=["GET"])
def form():
	return render_template("index.html")


@app.route("/run", methods=["POST"])
def run():
	allele = request.files["allele"]
	response = request.files["response"]
	prediction = request.files.get("prediction")

	run_id = datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:6]
	run_dir = UPLOAD_DIR / run_id
	run_dir.mkdir()

	allele_path = run_dir / allele.filename
	response_path = run_dir / response.filename
	allele.save(allele_path)
	response.save(response_path)

	out_dir = run_dir / "output"
	out_dir.mkdir()

	cmd = [
		sys.executable,
		str(SCRIPT_PATH),
		"-a", str(allele_path),
		"-r", str(response_path),
		"-o", str(out_dir)
	]

	if prediction and prediction.filename:
		pred_path = run_dir / prediction.filename
		prediction.save(pred_path)
		cmd += ["-p", str(pred_path)]

	response_cutoff = request.form.get("response_cutoff") or "1"
	cmd += ["-k", response_cutoff]

	if request.form.get("prediction_cutoff"):
		cmd += ["-c", request.form["prediction_cutoff"]]

	subprocess.run(cmd, capture_output=True, text=True)

	zip_path = shutil.make_archive(str(out_dir), "zip", root_dir=str(out_dir))
	Path(zip_path).replace(run_dir / "result.zip")

	# return render_template("ready.html", run_id=run_id)
	return render_template("running.html", run_id=run_id)


@app.route("/ready/<run_id>")
def ready(run_id):
    return render_template("ready.html", run_id=run_id)


@app.route("/ready/<run_id>", methods=["GET"])
def ready(run_id):
	return render_template("ready.html", run_id=run_id)


@app.route("/download/<run_id>", methods=["GET"])
def download(run_id):
	run_dir = UPLOAD_DIR / run_id
	zip_file = run_dir / "result.zip"

	response = send_file(
		zip_file,
		as_attachment=True,
		download_name=f"output_{run_id}.zip"
	)

	shutil.rmtree(run_dir, ignore_errors=True)
	return response


app.run(host="0.0.0.0", port=8000, debug=True)
