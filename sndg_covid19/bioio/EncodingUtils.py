

class EncodingUtils:

    DEFAULT_ENCODIG_LIST = tuple(["iso-8859-1","utf-8"])

    @staticmethod
    def get_encoding(file_path,encoding_list=DEFAULT_ENCODIG_LIST):
        try:
            with open(file_path,encoding=encoding_list[0]) as h:
                for x in h:
                    pass
            return encoding_list[0]
        except  UnicodeDecodeError:
            if len(encoding_list) > 1:
                return EncodingUtils.validate_endoding(file_path,encoding_list[1:])
            else:
                return None
