from aiosqlite import Connection
from argon2 import PasswordHasher
from fastapi import Request, Response

import dotenv
import json
import jwt
import os
import time

dotenv.load_dotenv()

EXPIRES = 24 * 60 * 60
JWT_SECRET = os.getenv("JWT_SECRET")
if JWT_SECRET is None or JWT_SECRET == "":
  raise ValueError("JWT_SECRET is not provided in the .env")

class User:
  id: int
  name: str
  username: str

def gen_error(code: int, msg: str) -> Response:
  return Response(
    status_code=code,
    content=json.dumps({"error": msg}),
    headers={"Content-Type": "application/json"}
  )

def is_protected_path(path: str) -> bool:
  splits = path.split("/")
  return splits[1] == "discover" or splits[1] == "results"

def is_username_valid(uname: str) -> bool:
    return all(c.isalnum() or c == "_" for c in uname)

async def validate_token(db: Connection, token: str | None) -> User | None:
  try:
    if token is None:
      print("token is None")
      return None
    
    if not token.startswith("Bearer "):
      print("token is not started with Bearer ")
      return None
    
    token = token.replace("Bearer ", "", 1)

    decoded = jwt.decode(token, JWT_SECRET, algorithms="HS256")
    sub = decoded["sub"]
    if type(sub) is not str:
      print("sub is not string")
      return None
    if not sub.isnumeric():
      print("sub is not string numeric")
      return None
    
    user = None
    cur = await db.execute("SELECT * FROM user WHERE id = ?", (int(sub),))
    rows = await cur.fetchall()
    for row in rows:
      user = User()
      user.name = row[1] # name field is 1st field
      user.username = row[2] # username field is 2nd field

    await cur.close()

    return user
  except Exception as e:
    print(e) # Debugging purpose
    return None

async def AuthMiddleware(request: Request, call_next):
  if request.method == "OPTIONS": # Biarin OPTIONS lewat
    return await call_next(request)

  if not hasattr(request.app.state, "db"):
    print("App state doesn't have db attribute")
    return gen_error(500, "Internal server error")

  db = request.app.state.db
  if not isinstance(db, Connection):
    print("Invalid db state")
    return gen_error(500, "Internal server error")

  token = request.headers.get("Authorization")
  validated = await validate_token(db, token)
  request.state.user = validated

  path = request.url.path

  if is_protected_path(path) and (validated is None):
    return gen_error(401, "Unauthorized")
  
  return await call_next(request)

async def LoginHandler(request: Request):
  if not hasattr(request.app.state, "db"):
    print("App state doesn't have db attribute")
    return gen_error(500, "Internal server error")

  db = request.app.state.db
  if not isinstance(db, Connection):
    print("Invalid db state")
    return gen_error(500, "Internal server error")

  try:
    body = await request.json()
  except:
    return gen_error(400, "Invalid request payload")

  uname = body.get("username", "")
  pw = body.get("password", "")
  if type(uname) is not str or type(pw) is not str:
    return gen_error(400, "username or password has invalid type (expected string)")
  
  uname = uname.strip()
  pw =  pw.strip()
  if len(uname) < 3:
    return gen_error(400, "Username length must greater than 3")
  if len(pw) < 8:
    return gen_error(400, "Password length must greater than 8")
  if not is_username_valid(uname):
    return gen_error(400, "Username must be a-z, A-Z, 0-9, and _")

  user_id, user_pass_hash = None, None
  cur = await db.execute("SELECT * FROM user WHERE username = ?", (uname,))
  rows = await cur.fetchall()
  for row in rows:
    user_id = row[0]
    user_pass_hash = row[3]
  
  if user_id is None or user_pass_hash is None:
    return gen_error(404, "User not found")

  hasher = PasswordHasher()
  verified = hasher.verify(user_pass_hash, pw)
  if not verified:
    return gen_error(401, "Wrong password")

  now = int(time.time())
  payload = {
    "sub": str(user_id),    # Subject, the user id, must be string, smh
    "nbf": now,             # Not before, ya, harus sekarang
    "iat": now,             # Kapan token dibuat, ya sekarang
    "exp": now + EXPIRES    # Kapan token basi
  }

  encoded = jwt.encode(payload, JWT_SECRET, algorithm="HS256")
  resp = { "token": encoded }

  return Response(
    status_code=200,
    content=json.dumps(resp),
    headers={"Content-Type":"application/json"}
  )

async def RegisterHandler(request: Request):
  if not hasattr(request.app.state, "db"):
    print("App state doesn't have db attribute")
    return gen_error(500, "Internal server error")

  db = request.app.state.db
  if not isinstance(db, Connection):
    print("Invalid db state")
    return gen_error(500, "Internal server error")

  try:
    body = await request.json()
  except:
    return gen_error(400, "Invalid request payload")

  name = body.get("name", "")
  uname = body.get("username", "")
  pw = body.get("password", "")
  if type(name) is not str or type(uname) is not str or type(pw) is not str:
    return gen_error(400, "username or password has invalid type (expected string)")
  
  name = name.strip()
  uname = uname.strip()
  pw =  pw.strip()
  if len(name) < 3:
    return gen_error(400, "Name length must greater than 3")
  if len(uname) < 3:
    return gen_error(400, "Username length must greater than 3")
  if len(pw) < 8:
    return gen_error(400, "Password length must greater than 8")
  if not is_username_valid(uname):
    return gen_error(400, "Username must be a-z, A-Z, 0-9, and _")
  
  user_exists = False
  cur = await db.execute("SELECT * FROM user WHERE username = ?", (uname,))
  rows = await cur.fetchall()
  for _ in rows:
    user_exists = True
  await cur.close()
  
  if user_exists:
    return gen_error(400, "Username is already taken")

  hasher = PasswordHasher()
  pass_hash = hasher.hash(pw)

  await db.execute("INSERT INTO user(name, username, password) VALUES(?, ?, ?)", (name, uname, pass_hash))
  await db.commit()

  return Response(
    status_code=201,
    content=json.dumps({"message": "Created."}),
    headers={"Content-Type":"application/json"}
  )