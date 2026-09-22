import pytest
from unittest.mock import patch, MagicMock
from src.download import ensure_data

@patch("src.download.REQUIRED_FILES", ["dummy_file.txt"])
@patch("src.download.DATA_DIR")
def test_ensure_data_already_present(mock_data_dir, capsys):
    # Mock that the file exists
    mock_file = MagicMock()
    mock_file.exists.return_value = True
    mock_data_dir.__truediv__.return_value = mock_file
    
    ensure_data()
    
    captured = capsys.readouterr()
    assert "[+] All required datasets are already present in data/" in captured.out

@patch("src.download.REQUIRED_FILES", ["dummy_missing_file.txt"])
@patch("src.download.DATA_DIR")
@patch("src.download.hf_hub_download")
@patch("src.download.login")
def test_ensure_data_missing(mock_login, mock_hf_download, mock_data_dir, capsys):
    # Mock that the file does NOT exist
    mock_file = MagicMock()
    mock_file.exists.return_value = False
    mock_data_dir.__truediv__.return_value = mock_file
    
    ensure_data()
    
    captured = capsys.readouterr()
    assert "[*] Missing files detected:" in captured.out
    assert "dummy_missing_file.txt" in captured.out
    
    mock_hf_download.assert_called_once()
