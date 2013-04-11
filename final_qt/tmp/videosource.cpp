#include "videosource.h"
#include <cstdlib>
#include <cstring>
#define VIDEO_WIDTH 640
#define VIDEO_HEIGHT 480

VideoSource::VideoSource()
{
    //preparation of the vlc command
    const char * const vlc_args[] = {
        "-I", "dummy", /* Don't use any interface */
        "--ignore-config", /* Don't use VLC's config */
        "--extraintf=logger", //log anything
        "--verbose=2", //be much more verbose then normal for debugging purpose
        "--plugin-path=C:\\vlc-0.9.9-win32\\plugins\\"
    };

    _isPlaying=false;

    //create a new libvlc instance
    //    _vlcinstance=libvlc_new(sizeof(vlc_args) / sizeof(vlc_args[0]), vlc_args,&_vlcexcep);  //tricky calculation of the char space used
    _vlcinstance=libvlc_new(sizeof(vlc_args) / sizeof(vlc_args[0]), vlc_args);
    //    raise (&_vlcexcep);

    // Create a media player playing environement
    //    _mp = libvlc_media_player_new (_vlcinstance, &_vlcexcep);
    _mp = libvlc_media_player_new (_vlcinstance);
    //    raise (&_vlcexcep);
}

//desctructor
VideoSource::~VideoSource()
{
    /* Stop playing */
    //    libvlc_media_player_stop (_mp, &_vlcexcep);
    libvlc_media_player_stop (_mp);

    /* Free the media_player */
    libvlc_media_player_release (_mp);

    libvlc_release (_vlcinstance);
    //    raise (&_vlcexcep);
}

void VideoSource::playFile(QString file)
{
    //the file has to be in one of the following formats /perhaps a little bit outdated)
    /*
    [file://]filename              Plain media file
    http://ip:port/file            HTTP URL
    ftp://ip:port/file             FTP URL
    mms://ip:port/file             MMS URL
    screen://                      Screen capture
    [dvd://][device][@raw_device]  DVD device
    [vcd://][device]               VCD device
    [cdda://][device]              Audio CD device
    udp:[[<source address>]@[<bind address>][:<bind port>]]
    */

    /* Create a new LibVLC media descriptor */
    //    _m = libvlc_media_new (_vlcinstance, file.toAscii(), &_vlcexcep);
    _m = libvlc_media_new_path (_vlcinstance, file.toAscii());
    //    raise(&_vlcexcep);

    // stuff for storing buffers
    struct ctx* context;
    context = (struct ctx*) malloc( sizeof(*context));
    context->pixels = (uchar *) malloc( (sizeof(*(context->pixels)) *VIDEO_WIDTH*VIDEO_HEIGHT) * 4);
    // mutex
    context->mutex = new QMutex();


    //    libvlc_media_player_set_media (_mp, _m, &_vlcexcep);
    libvlc_media_player_set_media (_mp, _m);
    //    raise(&_vlcexcep);

    // /!\ Please note /!\
    //
    // passing the widget to the lib shows vlc at which position it should show up
    // vlc automatically resizes the video to the given size of the widget
    // and it even resizes it, if the size changes at the playing

    /* Play */
    //    libvlc_media_player_play (_mp, &_vlcexcep );
    libvlc_media_player_play (_mp);
    //    raise(&_vlcexcep);

    _isPlaying=true;
}

void VideoSource::changeVolume(int newVolume)
{
    //    libvlc_exception_clear(&_vlcexcep);
    //    libvlc_audio_set_volume (_vlcinstance,newVolume , &_vlcexcep);
    //    libvlc_audio_set_volume (_vlcinstance,newVolume);
    libvlc_audio_set_volume(_mp, newVolume);
    //    raise(&_vlcexcep);
}

void VideoSource::changePosition(int newPosition)
{
    //    libvlc_exception_clear(&_vlcexcep);
    // It's possible that the vlc doesn't play anything
    // so check before
    //    libvlc_media_t *curMedia = libvlc_media_player_get_media (_mp, &_vlcexcep);
    libvlc_media_t *curMedia = libvlc_media_player_get_media (_mp);
    //    libvlc_exception_clear(&_vlcexcep);
    if (curMedia == NULL)
        return;

//    float pos=(float)(newPosition)/(float)POSITION_RESOLUTION;
    //    libvlc_media_player_set_position (_mp, pos, &_vlcexcep);
//    libvlc_media_player_set_position (_mp, pos);
    //    raise(&_vlcexcep);
}

void VideoSource::updateInterface()
{
    if(!_isPlaying)
        return;

    // It's possible that the vlc doesn't play anything
    // so check before
    //    libvlc_media_t *curMedia = libvlc_media_player_get_media (_mp, &_vlcexcep);
    libvlc_media_t *curMedia = libvlc_media_player_get_media (_mp);
    //    libvlc_exception_clear(&_vlcexcep);
    if (curMedia == NULL)
        return;

    //    float pos=libvlc_media_player_get_position (_mp, &_vlcexcep);
    float pos=libvlc_media_player_get_position (_mp);
//    int siderPos=(int)(pos*(float)(POSITION_RESOLUTION));
    //    int volume=libvlc_audio_get_volume (_vlcinstance,&_vlcexcep);
    //    int volume=libvlc_audio_get_volume (_vlcinstance);
    int volume = libvlc_audio_get_volume(_mp);
}


/*
void Player::raise(libvlc_exception_t * ex)
{
    if (libvlc_exception_raised (ex))
    {
         fprintf (stderr, "error: %s\n", libvlc_exception_get_message(ex));
         exit (-1);
    }
}
//*/
