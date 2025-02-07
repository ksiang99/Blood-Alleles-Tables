import os
import requests

def download_file(url, save_dir):
    try:
        # Create a session to persist certain parameters across requests
        with requests.Session() as session:
            response = session.get(url, stream=True)
            response.raise_for_status()  # Check for HTTP errors

            # Get the file name from the Content-Disposition header if available
            content_disposition = response.headers.get('Content-Disposition')
            if content_disposition:
                file_name = content_disposition.split('filename=')[-1].strip('"')
            else:
                # Fall back to the last part of the URL if no header is available
                file_name = os.path.join(save_dir, url.split("/")[-1] + ".pdf")  # Assuming it's a PDF

            file_path = os.path.join(save_dir, file_name)

            # Check if file already exists
            if os.path.exists(file_path):
                print(f"File already exists: {file_path}. Skipping download.")
                return

            # Write the content to a file
            with open(file_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)

            print(f"Downloaded: {file_path}")

            # Rename the files
            isbt_index = file_name.find('ISBT')
            if isbt_index >= 0:
                new_file_name = file_name[isbt_index + len('ISBT'):]
                new_file_name = new_file_name.lstrip('-')
                new_file_path = os.path.join(save_dir, new_file_name)
                os.rename(file_path, new_file_path)

    except Exception as e:
        print(f"Failed to download {url}. Reason: {e}")

def download_files_from_links(links, save_dir):
    os.makedirs(save_dir, exist_ok=True)
    print(f"Files will be saved to: {save_dir}")

    for link in links:
        print(f"Processing: {link}")
        download_file(link, save_dir)
    
    print("ISBT Blood Tables Downloaded")

def read_links_from_file(file_path):
    links = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                # Strip whitespace and ignore comments
                stripped_line = line.strip()
                if stripped_line and not stripped_line.startswith("#"):
                    links.append(stripped_line.split('#')[0].strip())  # Remove comments at the end of the line
        return links
    
    except Exception as e:
        print(f"Failed to read links from file: {e}")
        return []

if __name__ == "__main__":
    # Path to the text file containing download links
    links_file_path = "ISBT_links.txt"
    
    # Read download links from the text file
    download_links = read_links_from_file(links_file_path)

    save_directory = "ISBT Blood Alleles Table"
    download_files_from_links(download_links, save_directory)