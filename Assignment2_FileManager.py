import os
import subprocess

class FileManager:
    def __init__(self, remote, ml_dir, mc_dir):
        if not all(isinstance(arg, str) for arg in [remote, ml_dir, mc_dir]):
            raise TypeError("All inputs have to be strings.")
        self.remote = remote.rstrip("/")
        self.ml_dir = ml_dir.rstrip("/")
        self.mc_dir = mc_dir.rstrip("/")

    def convertCloudToLocal(self, filename):
        if not filename.startswith(self.remote):
            raise ValueError(f"Invalid cloud path: {filename}. This has to start with the Remote name followed by a colon")
        Sole_file_name = os.path.basename(filename)
        local_path = os.path.join(self.ml_dir, Sole_file_name)
        return local_path

    def convertLocalToCloud(self, filename):
        if not filename.startswith(self.ml_dir):
            raise ValueError(f"Invalid local path: {filename}. This has to start with {self.ml_dir}")
        Sole_file_name = os.path.basename(filename)
        cloud_path = f"{self.remote}/{self.mc_dir}/{Sole_file_name}"
        return cloud_path

    def uploadData(self, filename):
        cloud_path = self.convertLocalToCloud(filename)
        cloud_path_dir = os.path.dirname(cloud_path)
        print(f"Uploading {filename} to {cloud_path_dir}")
        
        try:
            subprocess.run(["rclone", "copy", filename, cloud_path_dir], check=True)
            print("Upload successful.")
        except subprocess.CalledProcessError:
            print("Upload failed.")

    def downloadData(self, filename):
        local_path = self.convertCloudToLocal(filename)
        local_dir = os.path.dirname(local_path)

        os.makedirs(local_dir, exist_ok=True)

        print(f"Downloading {filename} to {local_path}")
        
        try:
            subprocess.run(["rclone", "copy", filename, local_dir], check=True)
            print("Download successful.")
        except subprocess.CalledProcessError:
            print("Download failed.")


remote = "DropboxRemote:/"
ml_dir = "/Users/stephenphillips/Documents/Work/Python/"
mc_dir = "Stephen Phillips/Stephen Phillips/Python Class/"

fm = FileManager(remote, ml_dir, mc_dir)

local_file = "/Users/stephenphillips/Documents/Work/Python/HighScoreMedium.txt"
cloud_file = "DropboxRemote:/Stephen Phillips/Stephen Phillips/Python Class/TrialDoc.txt"

print("Cloud to Local:", fm.convertCloudToLocal(cloud_file))
print("Local to Cloud:", fm.convertLocalToCloud(local_file))


fm.uploadData(local_file)
fm.downloadData(cloud_file)
