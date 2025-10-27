"""
Download macroeconomic data from the Brazilian Central Bank API.
This script downloads various Brazilian macroeconomic time series from the BCB (Banco
Central do Brasil) API and combines them into a single CSV file for analysis.
Series Downloaded:
    - IPCA: Consumer price inflation, monthly % change (code 433)
    - USD/BRL: Exchange rate, PTAX buy rate (code 3698)
    - M2: Money supply in BRL millions (code 27789)
Output:
    - CSV file saved to: lectures/lecture8/brazil_macro_data.csv
    - Contains daily/monthly observations from 2000 to present
    - Index: Date (datetime column)
    - Columns: One column per series
Functions:
    download_bcb_series: Downloads a single time series from BCB API
Usage:
    Run the script directly to download and save the data:
    $ python download_macro_data_lecture8.py
Dependencies:
    - pandas: Data manipulation
    - requests: HTTP requests to BCB API
    - pyprojroot: Project root directory detection
Notes:
    - The script handles missing data and errors gracefully
    - Empty or failed downloads are skipped with warnings
    - All series are aligned by date in the final DataFrame
"""

import pandas as pd
import requests
from datetime import datetime
from pyprojroot import here

dir_to_save = here() / 'lectures' / 'lecture8'

def download_bcb_series(series_code, start_date='01/01/2000', end_date=None):
    """Download time series from Brazilian Central Bank API"""
    if end_date is None:
        end_date = datetime.now().strftime('%d/%m/%Y')
    
    url = f'https://api.bcb.gov.br/dados/serie/bcdata.sgs.{series_code}/dados'
    params = {
        'formato': 'json',
        'dataInicial': start_date,
        'dataFinal': end_date
    }
    
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        
        # Check if response has content
        if not response.text or response.text.strip() == '':
            print(f"  Warning: Empty response for series {series_code}")
            return None
            
        data = response.json()
        
        if not data:
            print(f"  Warning: No data returned for series {series_code}")
            return None
        
        df = pd.DataFrame(data)
        df['data'] = pd.to_datetime(df['data'], format='%d/%m/%Y')
        df['valor'] = pd.to_numeric(df['valor'], errors='coerce')
        df = df.set_index('data')
        
        return df['valor']
        
    except Exception as e:
        print(f"  Error downloading series {series_code}: {e}")
        return None

# Series codes with descriptions
series_info = {
    'ipca': (433, 'IPCA - Consumer price inflation (monthly % change)'),
    'usd_brl': (3698, 'USD/BRL - Exchange rate (BRL per USD, daily)'),
    'm2': (27789, 'M2 - Money supply (BRL millions)')
}

# Extract codes for downloading
series = {name: code for name, (code, desc) in series_info.items()}

# Download all series
print("Downloading data from BCB...")
data = {}
for name, code in series.items():
    print(f"Downloading {name} (code {code})...", end=' ')
    series_data = download_bcb_series(code)
    if series_data is not None:
        data[name] = series_data
        print("✓")
    else:
        print("✗ (skipped)")

# Combine into single DataFrame
if data:
    df = pd.DataFrame(data)
    df.index.name = 'Date'
    
    # Save to CSV
    df.to_csv(dir_to_save / 'brazil_macro_data.csv')
    print(f"\nSaved to {dir_to_save / 'brazil_macro_data.csv'}")
    print(f"Observations: {len(df)}")
    print(f"Date range: {df.index.min()} to {df.index.max()}")
    print(f"Series downloaded: {list(data.keys())}")
else:
    print("\nNo data was successfully downloaded!")