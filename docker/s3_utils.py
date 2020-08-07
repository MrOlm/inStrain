import io
import os
import subprocess
import shlex
import boto3

s3 = boto3.resource('s3')

def download_folder(s3_path, directory_to_download):
    """
    Downloads a folder from s3
    :param s3_path: s3 folder path
    :param directory_to_download: path to download the directory to
    :return: directory that was downloaded
    """
    cmd = 'aws s3 cp --recursive %s %s' % (s3_path, directory_to_download)
    print(cmd)

    subprocess.check_call(shlex.split(cmd))

    return directory_to_download


def download_file(s3_path, directory_to_download):
    """
    Downloads an object from s3 to a local path
    :param s3_path: s3 object path
    :param directory_to_download: directory to download to
    :return: local file path of the object
    """
    bucket = s3_path.split('/')[2]
    key = '/'.join(s3_path.split('/')[3:])

    object_name = key.split('/')[-1]

    local_file_name = os.path.join(directory_to_download, object_name)

    s3.Object(bucket, key).download_file(local_file_name)

    return local_file_name


def upload_folder(s3_path, local_folder_path, sse=False):
    """
    Uploads a local folder to S3
    :param s3_path: s3 path to upload folder to
    :param local_folder_path: local folder path
    :param sse: boolean whether to enable server-side encryption
    """
    cmd = 'aws s3 cp --recursive %s %s' % (local_folder_path, s3_path)

    if sse:
        cmd += ' --sse'

    subprocess.check_call(shlex.split(cmd))


def upload_file(s3_path, local_path):
    """
    Uploads a local file to s3 with server side encryption enabled
    :param s3_path: s3 object path
    :param local_path: local file path
    :return: response from the upload file
    """
    bucket = s3_path.split('/')[2]
    key = '/'.join(s3_path.split('/')[3:])

    response = s3.Object(bucket, key).upload_file(local_path, ExtraArgs=dict(ServerSideEncryption='AES256'))

    return response


def get_matching_s3_objects(bucket, prefix="", suffix=""):
    """
    Generate objects in an S3 bucket.

    :param bucket: Name of the S3 bucket.
    :param prefix: Only fetch objects whose key starts with
        this prefix (optional).
    :param suffix: Only fetch objects whose keys end with
        this suffix (optional).
    """
    s3 = boto3.client('s3')
    paginator = s3.get_paginator("list_objects_v2")

    kwargs = {'Bucket': bucket}

    # We can pass the prefix directly to the S3 API.  If the user has passed
    # a tuple or list of prefixes, we go through them one by one.
    if isinstance(prefix, str):
        prefixes = (prefix, )
    else:
        prefixes = prefix

    for key_prefix in prefixes:
        kwargs["Prefix"] = key_prefix

        for page in paginator.paginate(**kwargs):
            try:
                contents = page["Contents"]
            except KeyError:
                return

            for obj in contents:
                key = obj["Key"]
                if key.endswith(suffix):
                    yield obj

def get_matching_s3_keys(bucket, prefix="", suffix=""):
    """
    Generate the keys in an S3 bucket.

    :param bucket: Name of the S3 bucket.
    :param prefix: Only fetch keys that start with this prefix (optional).
    :param suffix: Only fetch keys that end with this suffix (optional).
    """
    for obj in get_matching_s3_objects(bucket, prefix, suffix):
        yield obj["Key"]

def check_s3_file(floc):
    '''
    Return True if exists and False if it does not
    '''
    bucket = floc.split('/')[2]
    prefix = '/'.join(floc.split('/')[3:])

    found = False
    for key in get_matching_s3_keys(bucket, prefix):
        if prefix in key:
            found = True
    return found

def store_s3_file(bucket, location, binary_string):
    s3 = boto3.resource('s3')
    object = s3.Object(bucket, location)
    object.put(Body=binary_string)


def load_coverage_report(s3_bucket, s3_key, sep='\t', names=None):
    '''
    https://towardsdatascience.com/web-scraping-html-tables-with-python-c9baba21059
    '''
    # Load the data from s3
    client = boto3.client("s3")
    obj = client.get_object(Bucket=s3_bucket, Key=s3_key)
    df = pd.read_csv(io.BytesIO(obj['Body'].read()), sep=sep, names=names)

    return df

def glob_s3(s3_path):
    bucket = s3_path.split('/')[2]
    key = '/'.join(s3_path.split('/')[3:])

    prefix = key.split('*')[0]
    suffix = key.split('*')[1]

    return ['s3://' + bucket + '/' + x for x in get_matching_s3_keys(bucket, prefix, suffix)]
