import os
import requests
from dotenv import load_dotenv
import json

# Load environment variables
load_dotenv()

def test_nvidia_genmol_api():
    # Get API credentials
    api_key = os.getenv('NVIDIA_GENMOL_API_KEY')
    endpoint = os.getenv('NVIDIA_GENMOL_ENDPOINT')

    print(f"API Endpoint: {endpoint}")
    print(f"API Key (first 4 chars): {api_key[:4]}...")

    # Prepare headers and payload
    headers = {
        'Authorization': f'Bearer {api_key}',
        'Content-Type': 'application/json'
    }

    payload = {
        'smiles': 'C124CN3C1.S3(=O)(=O)CC.C4C#N.[*{20-20}]',
        'num_molecules': '30',
        'temperature': '1',
        'noise': '0',
        'step_size': '1',
        'scoring': 'QED'
    }

    try:
        # Send POST request
        response = requests.post(
            endpoint, 
            headers=headers, 
            data=json.dumps(payload)
        )

        # Print detailed response information
        print("\nResponse Details:")
        print(f"Status Code: {response.status_code}")
        print("Headers:")
        for key, value in response.headers.items():
            print(f"  {key}: {value}")
        
        print("\nResponse Content:")
        print(response.text)

        # Raise an exception for HTTP errors
        response.raise_for_status()

    except requests.exceptions.RequestException as e:
        print(f"\nRequest Error: {e}")
        print(f"Full Error Details: {response.text}")

if __name__ == '__main__':
    test_nvidia_genmol_api()
