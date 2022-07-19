def prefix_function(s):
    s_len = len(s)
    pi = [0] * s_len
    for i in range(1,s_len):
        j = pi[i-1]
        while j > 0 and s[i] != s[j]:
                j -= 1
        if s[i] == s[j]:
            pi[i] = j + 1
    return pi

class node:
    def __init__(self, char):
        self.char = char
        self.children = dict()
        self.count = 0
        self.parent = None
        self.is_end = False
        self.fail = None
        self.tag = None

class Automaton:
    def __init__(self) -> None:
        self.root = node(None)
        self.root.fail = self.root
        self.root.parent = node(None)
        self.root.parent.fail = self.root

    def insert(self, word, tag=None):
        root = self.root
        for idx, char in enumerate(word):
            if char in root.children:
                root = root.children[char]
                continue
            node_char = node(char)
            node_char.parent = root 
            root.children[char] = node_char
            root.count +=1
            if idx == 0:
                node_char.fail = self.root
            root = node_char
        root.tag = tag
        root.is_end = True
    
    def build(self):
        bfs_queue = []
        root_queue = [self.root.children[child] for child in self.root.children]
        '''
        for node in root_queue:
            for char in node.children:
                bfs_queue.append(node.children[char])
        root_queue = bfs_queue
        '''
        while root_queue:
            while root_queue:
                root = root_queue.pop(0)
                for child in root.children:
                    fail = root.fail
                    bfs_queue.append(root.children[child])
                    if child in fail.children:
                        root.children[child].fail = fail.children[child]
                        continue
                    while fail != self.root and (child not in fail.children):
                        fail = fail.fail
                    if child in fail.children:
                        root.children[child].fail = fail.children[child]
                    else:
                        root.children[child].fail = self.root
            root_queue = bfs_queue
            bfs_queue = []

    def print(self):
        bfs_queue = []
        root_queue = [self.root]
        while root_queue:
            while root_queue:
                root = root_queue.pop(0)
                parent = root.parent
                #print(root, root.char, root.parent, root.fail)
                if root.is_end:
                    print(root.word)
                for child in root.children:
                    bfs_queue.append(root.children[child])
                
            root_queue = bfs_queue
            bfs_queue = []

    def query(self, s):
        u = 0
        match = 0
        root = self.root

        for i in range(len(s)):
            root = self.root
            pos = i
            fail = False
            while pos < len(s):
                if s[pos] in root.children:
                    root = root.children[s[pos]]
                    pos += 1
                    if root.is_end and not fail:
                        match += 1
                        print(root.tag)
                else:
                    fail = True
                    root = root.fail
                    if root == self.root:
                        fail = False
                        break
        
                