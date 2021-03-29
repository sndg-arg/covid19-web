class JobValidationError(Exception):
    def __init__(self, status_text,errors):
        self.status_text = status_text
        self.errors = errors
        super().__init__(self.status_text)
