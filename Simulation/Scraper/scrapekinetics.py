import requests
from lxml import html
import csv
    
def search_database(order):
    search_url = 'http://kinetics.nist.gov/kinetics/Search.jsp'
    
    params = {
        'doc': 'SearchForm',
        'type': 'java',
        'Units':'',
        'database':'kinetics',
        'numberOfFields':5,
        'boolean1': '',
        'lp1': '+',
        'field1': 'rxn_order',
        'relate1': '=',
        'text1': str(order),
        'rp1': '+',
        'boolean2': 'and',
        'lp2': '+',
        'field2': 'reactants',
        'relate2': '=',
        'text2': '',
        'rp2': '+',
        'boolean3': 'and',
        'lp3': '+',
        'field3': 'products',
        'relate3': '=',
        'text3': '',
        'rp3': '+',
        'boolean4': 'and',
        'lp4': '+',
        'field4': 'products',
        'relate4': '=',
        'text4': '',
        'rp4': '+',
        'boolean5': 'and',
        'lp5': '+',
        'field5': 'kinetics.squib',
        'relate5': '~*',
        'text5': '',
        '&rp5': '+',
        'category': 0
        }
    
    page = requests.post(search_url, data=params)
    return page

def parse_reaction_page(page, csvwriter):
    tree = html.fromstring(page.content)
    #<div align="center"><font size="+3"><b>C<sub>2</sub>H<sub>5</sub>OCH=CHNH<sub>2</sub> → <a href="http://webbook.nist.gov/cgi/cbook.cgi?ID=74851&amp;Units=SI" target=_blank onMouseOver="showPicture('http://webbook.nist.gov/cgi/cbook.cgi?Struct=74851',event);" onMouseOut="Hide('iAlt');">C<sub>2</sub>H<sub>4</sub></a> + <a href="http://webbook.nist.gov/cgi/cbook.cgi?ID=74895&amp;Units=SI" target=_blank onMouseOver="showPicture('http://webbook.nist.gov/cgi/cbook.cgi?Struct=74895',event);" onMouseOut="Hide('iAlt');">CH<sub>3</sub>NH<sub>2</sub></a> + <a href="http://webbook.nist.gov/cgi/cbook.cgi?ID=630080&amp;Units=SI" target=_blank onMouseOver="showPicture('http://webbook.nist.gov/cgi/cbook.cgi?Struct=630080',event);" onMouseOut="Hide('iAlt');">CO</a></b></font></div>
    reaction_pieces = tree.xpath('''//div[@align="center"]/font/b/text()
        |//div[@align="center"]/font/b/sub/text()
        |//div[@align="center"]/font/b/a/text()
        |//div[@align="center"]/font/b/a/sub/text()''')
    reaction_string = ''.join(reaction_pieces)
    reaction_string = reaction_string.replace('&middot', '.')
    reaction_string = reaction_string.replace('·', '.')
    
    if 'Products' in reaction_string:
        return #skip reactions without listed products
    try:
        (reactants, products) = reaction_string.split('→')
    except ValueError as e:
        print('Failed to split: %s' % reaction_string)
        return
    #Split each side into lists of molecules
    reactants = reactants.split('+')
    products = products.split('+')
    
    #Remove whitespace
    reactants = [x.strip() for x in reactants]
    products = [x.strip() for x in products]
    
    print('%s -> %s'% (str(reactants), str(products)))

    data_rows = tree.xpath('//table[contains(*, "Temp")]/tr[position() > 2]')
    #print(data_rows)
    success = False
    for dataset in data_rows:
        #print(html.tostring(dataset))
        try:
            #Temp
            TempRange = dataset.xpath('./td[5]/text()')
            if not TempRange:
                continue
            TempRange = TempRange[0].split('-')
            TempRange = [float(t) for t in TempRange] #Convert to float
            if len(TempRange) == 1:
                TempRange.append(TempRange[0])
            (LowTemp, HighTemp) = min(TempRange), max(TempRange)

            #A
            A = float(dataset.xpath('./td[7]/text()')[0])

            #n - optional
            n = dataset.xpath('./td[9]/text()')[0]

            #Ea
            Ea = float(dataset.xpath('./td[11]/text()')[0])
            
            #Order
            Order = int(dataset.xpath('./td[15]/text()')[0])

        except ValueError as e:
            print(e)
            continue #Empty columns, try next row
        else:
            success = True
            break #Success, don't try next row
    #End for loop    
    if not success:
        return

    ascii_reaction_string = '+'.join(reactants) + '->' + '+'.join(products)

    Rall = ' | '.join(reactants)
    r_max = 2;
    while len(reactants) < r_max:
        reactants.append('')
        
    Pall = ' | '.join(products)
    p_max = 4;
    while len(products) < p_max:
        products.append('')
    
    print('LowTemp=%0.1f HighTemp=%0.1f A=%G n=%s Ea=%0.1f Order=%d' % (LowTemp, HighTemp, A, n, Ea, Order))
    try:
        csvwriter.writerow( [ascii_reaction_string, Order, Rall, reactants[0], reactants[1], Pall, products[0], products[1], products[2], products[3], LowTemp, HighTemp, A, n, Ea] )
    except UnicodeEncodeError as e:
        pass #ignore rows which contain non-ascii characters
    
def parse_result_page(page, csvwriter):
    tree = html.fromstring(page.content)

    #Build list of reaction page links:
    #<tr bgcolor="#e8e8e8"><td><a href="/kinetics/ReactionSearch?r0=1154406344&r1=0&r2=0&r3=0&r4=0&p0=74851&p1=74895&p2=630080&p3=0&p4=0&expandResults=true&">1 record matched</a></td><td width=10>&nbsp;</td><td>C<sub>2</sub>H<sub>5</sub>OCH=CHNH<sub>2</sub> → C<sub>2</sub>H<sub>4</sub> + CH<sub>3</sub>NH<sub>2</sub> + CO</td></tr>
    links = tree.xpath('//table/tr/td[1]/a/@href')
    for l in links:
        reaction_link = 'http://kinetics.nist.gov' + l;
        reaction_page = requests.get(reaction_link)
        parse_reaction_page(reaction_page, csvwriter)
    
def scrape():
    #open output CSV file
    with open('reactions2.csv', 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
        csvwriter.writerow(['Reaction', 'Order', 'Rall', 'R1', 'R2', 'Pall', 'P1', 'P2', 'P3', 'P4', 'Tmin', 'Tmax', 'A', 'n', 'Ea'])
    
        first_order_results = search_database(order=1)
        parse_result_page(first_order_results, csvwriter)
        second_order_results = search_database(order=2)
        parse_result_page(second_order_results, csvwriter)


def main():
    scrape()
    
if __name__ == "__main__":
    main()

    
