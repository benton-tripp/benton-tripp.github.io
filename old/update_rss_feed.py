import os
from xml.etree.ElementTree import Element, SubElement, tostring, ElementTree
from xml.dom import minidom

def prettify(elem):
    rough_string = tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def generate_rss():
    rss = Element('rss', version='2.0')
    channel = SubElement(rss, 'channel')
    title = SubElement(channel, 'title')
    title.text = "Benton's Portfolio/Blog"
    
    for filename in os.listdir('posts'):
        if filename.endswith('.html'):
            item = SubElement(channel, 'item')
            item_title = SubElement(item, 'title')
            item_title.text = filename.replace('.html', '')
            link = SubElement(item, 'link')
            link.text = f'https://benton-tripp.github.io/posts/{filename}'
    
    tree = ElementTree(rss)
    with open("feed.xml", "wb") as fh:
        tree.write(fh, encoding='utf-8', xml_declaration=True)

    print('RSS feed generated as feed.xml')

if __name__ == "__main__":
    generate_rss()
