"""Render the Markdown paper and supplement with Pandoc and their embedded figures."""
from pathlib import Path
from tempfile import TemporaryDirectory
import re
import shutil
import subprocess
import sys

import cairosvg
from docx import Document
from docx.shared import Inches, Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH

assets = Path(__file__).resolve().parent
output = Path(sys.argv[1]).resolve()
output.mkdir(parents=True, exist_ok=True)
sources = [(assets.parent/'FCP_NEW_PHYTOLOGIST_INTEGRATED_POSITIVE_RESULTS_20261024.md', 'FCP_Manuscript', 3),
           (assets/'Supporting_Information.md', 'FCP_Supporting_Information', 2)]
for source, name, count in sources:
    with TemporaryDirectory() as folder:
        tmp = Path(folder)
        def prepare_image(match):
            alt, link = match.groups()
            path = (source.parent/link).resolve()
            if path.suffix == '.svg':
                png = tmp/(path.stem+'.png')
                cairosvg.svg2png(url=str(path), write_to=str(png), output_width=1500)
                path = png
            return f'![{alt}]({path.as_posix()}){{width=6in}}'
        body = re.sub(r'!\[([^\]]*)\]\(([^)]+)\)', prepare_image, source.read_text(encoding='utf8'))
        # The Word and PDF supplement travels next to the manuscript.
        body = body.replace('(climate_colour/Supporting_Information.md)', '(FCP_Supporting_Information.pdf)')
        md = tmp/'input.md'
        md.write_text(body, encoding='utf8')
        file = output/(name+'.docx')
        subprocess.run(['pandoc',str(md),'-f','markdown-implicit_figures','-o',str(file)],check=True)
        doc=Document(file)
        doc.core_properties.title=source.read_text(encoding='utf8').splitlines()[0].lstrip('# ')
        for section in doc.sections:
            section.top_margin=section.bottom_margin=Inches(.8)
            section.left_margin=section.right_margin=Inches(.95)
        for index, paragraph in enumerate(doc.paragraphs):
            if paragraph._p.xpath('.//w:drawing'):
                paragraph.alignment=WD_ALIGN_PARAGRAPH.CENTER
                paragraph.paragraph_format.page_break_before=True
                paragraph.paragraph_format.keep_with_next=True
            if paragraph.text.startswith(('Figure ','Figure S')):
                paragraph.paragraph_format.keep_together=True
            if paragraph.text=='Figures':
                # Use the first caption as the figure section opener.
                paragraph._p.getparent().remove(paragraph._p)
        doc.styles['Normal'].font.size=Pt(11)
        doc.save(file)
        assert len(Document(file).inline_shapes)==count
for filename in ['Table_S1_photographs.csv','Table_S2_species_coefficients.csv.gz']:
    shutil.copyfile(assets/filename,output/filename)
