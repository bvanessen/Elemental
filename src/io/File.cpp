/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {

const char* QtImageFormat( FileFormat format )
{
    switch( format )
    {
    case FileFormat::BMP:  return "FileFormat::BMP";  break;
    case FileFormat::JPG:  return "JPG";  break;
    case FileFormat::JPEG: return "JPEG"; break;
    case FileFormat::PNG:  return "PNG";  break;
    case FileFormat::PPM:  return "PPM";  break;
    case FileFormat::XBM:  return "XBM";  break;
    case FileFormat::XPM:  return "XPM";  break;
    default: LogicError("Invalid image format"); return "N/A"; break;
    }
}

string FileExtension( FileFormat format )
{
    switch( format )
    {
    case FileFormat::ASCII:            return "txt";  break;
    case FileFormat::ASCII_MATLAB:     return "m";    break;
    case FileFormat::BINARY:           return "bin";  break;
    case FileFormat::BINARY_FLAT:      return "dat";  break;
    case FileFormat::BMP:              return "bmp";  break;
    case FileFormat::JPG:              return "jpg";  break;
    case FileFormat::JPEG:             return "jpeg"; break;
    case FileFormat::MATRIX_MARKET:    return "mm";   break;
    case FileFormat::PNG:              return "png";  break;
    case FileFormat::PPM:              return "ppm";  break;
    case FileFormat::XBM:              return "xbm";  break;
    case FileFormat::XPM:              return "xpm";  break;
    default: LogicError("Format not found"); return "N/A"; break;
    }
}

FileFormat FormatFromExtension( const string ext )
{
    bool foundFormat = false;
    FileFormat format = FileFormat::BINARY;
    for( int j=1; j<FileFormat_MAX; ++j )
    {
        format = static_cast<FileFormat>(j);
        if( FileExtension(format) == ext )
        {
            foundFormat = true;
            break;
        }
    }
    if( !foundFormat )
        RuntimeError("Did not detect file format");
    return format;
}

FileFormat DetectFormat( const string filename )
{
    const string ext = filename.substr(filename.find_last_of(".")+1);
    return FormatFromExtension( ext );
}

ifstream::pos_type FileSize( ifstream& file )
{
    auto pos = file.tellg();
    file.seekg( 0, ifstream::end );
    auto numBytes = file.tellg();
    file.seekg( pos );
    return numBytes;
}

} // namespace El
