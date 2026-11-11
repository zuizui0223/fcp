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
from docx.oxml import OxmlElement

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
                paragraph.paragraph_format.page_break_before=(name=='FCP_Manuscript')
                paragraph.paragraph_format.keep_with_next=True
            if paragraph.text.startswith(('Figure ','Figure S')):
                paragraph.paragraph_format.keep_together=True
            if paragraph.text=='Figures':
                # Use the first caption as the figure section opener.
                paragraph._p.getparent().remove(paragraph._p)
        for table in doc.tables:
            widths = [0.45, 2.1, 1.5, 0.65, 0.65, 0.65] if len(table.columns)==6 else [2.6, 0.8, 2.6]
            table.autofit=False
            for col, width in zip(table.columns,widths):
                col.width=Inches(width)
                for cell in col.cells:
                    cell.width=Inches(width)
            for row in table.rows:
                row._tr.get_or_add_trPr().append(OxmlElement('w:cantSplit'))
                for cell in row.cells:
                    for para in cell.paragraphs:
                        para.paragraph_format.space_after=Pt(2)
                        for run in para.runs:
                            run.font.size=Pt(9)
        doc.styles['Normal'].font.size=Pt(11)
        doc.save(file)
        assert len(Document(file).inline_shapes)==count
for filename in ['Table_S1_photographs.csv','Table_S2_species_coefficients.csv.gz']:
    shutil.copyfile(assets/filename,output/filename)
