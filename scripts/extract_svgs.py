#!/usr/bin/env python3
import os
import sys
import base64
from bs4 import BeautifulSoup

def extract_svgs(html_file, target_dir):
    with open(html_file, 'r', encoding='utf-8') as f:
        soup = BeautifulSoup(f, 'html.parser')
    
    # Find the container with the expected id.
    container = soup.select_one('#visualize-suitable-sampling-regions')
    if container is None:
        print("Container with id 'visualize-suitable-sampling-regions' not found.")
        return

    # Look for div elements whose id ends with '_bioclim_plots'
    plot_divs = container.find_all('div', id=lambda x: x and x.endswith('_bioclim_plots'))
    if not plot_divs:
        print("No plot divs found within the container.")
        return

    for div in plot_divs:
        div_id = div.get('id')
        # Expecting an id of the form "<state>_<bird>_bioclim_plots"
        parts = div_id.split('_')
        if len(parts) < 3:
            print(f"Skipping div with unexpected id format: {div_id}")
            continue

        state = parts[0]
        bird = parts[1]
        imgs = div.find_all('img')
        for i, img in enumerate(imgs, start=1):
            src = img.get('src', '')
            if src.startswith("data:"):
                try:
                    header, b64data = src.split(',', 1)
                except ValueError:
                    print(f"Skipping image in div {div_id} with invalid data URI.")
                    continue
                try:
                    svg_data = base64.b64decode(b64data)
                except Exception as e:
                    print(f"Error decoding base64 data in div {div_id}: {e}")
                    continue

                file_name = f"sdm-6-{state}-{bird}-{i}.svg"
                file_path = os.path.join(target_dir, file_name)
                with open(file_path, 'wb') as out_file:
                    out_file.write(svg_data)
                print(f"Saved {file_path}")
            else:
                print(f"Skipping non-data image in div {div_id}")

if __name__ == "__main__":
    # if len(sys.argv) < 3:
    #     print("Usage: python extract_svgs.py input.html target_directory")
    #     sys.exit(1)

    html_file = "../old/posts/2023-10-18-sdm-benchmark-study-part-6-resampling-pseudoabsence-points.html "# sys.argv[1]
    target_dir = "../assets/plots" # sys.argv[2]
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    extract_svgs(html_file, target_dir)
