import json
import os

from jinja2 import Environment, FileSystemLoader

# Load JSON file
with open("data/pango-consensus-sequences_summary.json", "r") as f:
    data = json.load(f)

# Create output directory if it doesn't exist
if not os.path.exists("output"):
    os.makedirs("output")

# Set up Jinja2 environment
env = Environment(loader=FileSystemLoader("."))
template = env.get_template("web/lineage_template.html")

# Loop through each lineage and generate HTML page
for lineage in data.keys():
    # Create lineage directory if it doesn't exist
    # Render Jinja2 template with lineage data
    output = template.render(lineage=data[lineage])
    # Write HTML file to lineage directory
    with open(f"output/{lineage}.html", "w") as f:
        f.write(output)
